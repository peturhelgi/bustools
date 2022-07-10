#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <vector>
#include <array>
#include <stdint.h>

#include <assert.h>
#include "Common.hpp"
#include "BUSData.h"
#include "bustools_compress.h"
#include "bustools_decompress.h"

/**
 * @brief Decode a single fibonacci number from a buffer
 * @pre buf contains at least n_buf elements.
 * @pre 0 <= i_fibo < fibo64.size()
 * @pre 0 <= bitpos < 64
 * @pre 0 <= bufpos < n_buf
 *
 * @param buf The buf to fibo-decode from.
 * @param n_buf maximum number of elements in the buffer.
 * @param i_fibo a counter of previously seen fibonacci numbers.
 * @param bitpos a bit offset from the left within buf[bufpos].
 * @param bufpos offset from the beginning of buf.
 * 
 * @return uint64_t num; the fibonacci decoded number. If buf is exhausted before the termination bit of
 * a fibonacci number, num represents a partially decoded number.
 */
uint64_t fiboDecodeSingle(uint64_t const *const buf, const size_t n_buf, uint32_t &i_fibo, uint32_t &bitpos, uint32_t &bufpos)
{

	size_t SIZE = sizeof(uint64_t) * 8;
	size_t buf_offset = bufpos;
	size_t bit_offset = SIZE - (bitpos % SIZE) - 1;
	uint64_t num{0};

	int bit = (buf[buf_offset] >> bit_offset) & 1;
	int last_bit = 0;

	++bitpos;
	bool print =0;
	bool print_fibo = print && 0;
	bool print_bits = print && !print_fibo;

	while (last_bit + bit < 2 && buf_offset < n_buf)
	{
		last_bit = bit;
		num += bit * (fibo64[i_fibo]);

		buf_offset = bufpos + bitpos / SIZE;
		bit_offset = SIZE - (bitpos % SIZE) - 1;

		bit = (buf[buf_offset] >> bit_offset) & 1;

		++bitpos;
		++i_fibo;
	}
	if (last_bit + bit == 2)
	{
		i_fibo = 0;
	}

	bufpos += (bitpos / SIZE);
	bitpos %= SIZE;
	return num;
}

/**
/**
 * @brief 
 * 
 * @param buf 
 * @param max_elem 
 * @param n_ints 
 * @param b_bits 
 * @param primary 
 * @param bit_pos 
 * @param buf_offset 
 * @return size_t 
 */
size_t PfdParsePrimaryBlock(uint64_t *buf, const size_t max_elem, const int n_ints, const uint32_t b_bits, int32_t *primary, uint32_t &bit_pos, size_t buf_offset)
{

	int i = 0;
	while (i < n_ints && buf_offset < max_elem)
	{
		// I: 0 <= bit_pos < 64

		uint32_t n_partial = std::min(b_bits, 64 - bit_pos);
		uint32_t bits_rem = b_bits - n_partial;

		uint32_t shift = (64 - n_partial - bit_pos);
		uint64_t mask = (1 << n_partial) - 1;

		int32_t num = ((buf[buf_offset] >> shift) & mask) << bits_rem;

		bit_pos = (bit_pos + n_partial) % 64;
		// increment buf_offset if bit_pos == 0, i.e. the next element in `buf` is reached.
		buf_offset += (!bit_pos);
		
		if (bits_rem)
		{
			shift = 64 - bits_rem;
			num |= buf[buf_offset] >> shift;
			bit_pos += bits_rem;
		}
		primary[i++] = num;
	}
	// std::cout << "i: " << i << '\n';
	return buf_offset;
}



/**
 * @brief Decompress barcodes using fibonacci-runlength(0)-delta decoding.
 * 
 * @param BUF The char array occupied by at least `row_count` encoded numbers, encoded using delta-runlenght(0)-fibonacci.
 * @param rows The output BUSData array.
 * @param row_count The number of BUSData elements in `rows`.
 * @param buf_size The size of `BUF` in bytes.
 */
void decompress_barcode(char *BUF, BUSData *rows, const size_t row_count, const size_t buf_size) {
	uint64_t *BUF64 = (uint64_t *)BUF;
	size_t buf64_size = (buf_size - 1) / (sizeof(uint64_t)) + 1;

	uint32_t bitpos{0},
		i_fibo{0},
		buf_offset{0};
	size_t row_index = 0;
	uint64_t diff = 0,
			 barcode = 0,
			 runlen = 0;

	uint64_t RLE_VAL{0ULL};
	while (row_index < row_count)
	{
		diff = fiboDecodeSingle(BUF64, buf64_size, i_fibo, bitpos, buf_offset) - 1;
		
		if (diff == RLE_VAL)
		{
			// Runlength decoding
			runlen = fiboDecodeSingle(BUF64, buf64_size, i_fibo, bitpos, buf_offset);
			for (int i = 0; i < runlen; ++i)
			{
				rows[row_index].barcode = barcode;
				++row_index;
			}
		}
		else
		{
			// Delta decoding
			barcode += diff;
			rows[row_index].barcode = barcode;
			++row_index;
		}
	}
}

/**
 * @brief Decompress UMIS using fibonacci-runlength-periodic_delta decoding.
 * @pre rows[0]..rows[row_count - 1] have its barcodes set from `decompress_barcodes`.
 * @pre rows have at least `row_count` elements.
 * 
 * @param BUF The encoded UMIs to decode.
 * @param rows The output BUSData array whose UMI members are set in this function (up to `row_count`-1).
 * @param row_count The number of elements in `rows` to decode.
 * @param buf_size The size of the input array `BUF`.
 */
void decompress_lossless_umi(char *BUF, BUSData *rows, const size_t row_count, const size_t buf_size) {
	uint64_t *BUF64 = (uint64_t *)BUF;
	size_t buf64_size = (buf_size - 1) / (sizeof(uint64_t)) + 1;

	uint32_t bitpos{0},
		i_fibo{0},
		buf_offset{0};
	size_t row_index = 0;
	uint64_t diff = 0,
			 last_barcode = rows[0].barcode + 1,
			 umi = 0,
			 barcode,
			 runlen = 0;

	uint64_t RLE_VAL{0ULL};

	while(row_index < row_count){
		diff = fiboDecodeSingle(BUF64, buf64_size, i_fibo, bitpos, buf_offset) - 1;
		barcode = rows[row_index].barcode;
		if (barcode != last_barcode)
		{
			umi = 0;
		}

		if (diff == RLE_VAL)
		{
			// diff is runlen encoded, next values are identical.
			runlen = fiboDecodeSingle(BUF64, buf64_size, i_fibo, bitpos, buf_offset);
			for (int i = 0; i < runlen; ++i)
			{
				rows[row_index].UMI = umi - 1;
				++row_index;
			}
		}
		else{
			// Current element is a delta
			umi += diff;
			rows[row_index].UMI = umi - 1;
			++row_index;
		}
		last_barcode = barcode;
	}
}
/**
 * @brief Decompress UMIs which have been compressed in a lossy manner.
 * @note Not yet implemented.
 * @pre rows[0]..rows[row_count - 1] have its barcodes set from `decompress_barcodes`.
 * @pre rows have at least `row_count` elements.
 * 
 * @param BUF The encoded UMIs to decode
 * @param rows The output BUSData array whose UMI members are set in this function (up to `row_count`-1).
 * @param row_count The number of elements in `rows` to decode.
 * @param buf_size The size of the input array `BUF`.
 */
void decompress_lossy_umi(char *BUF, BUSData *rows, const size_t row_count, const size_t buf_size) {}

/**
 * @brief Decompress ECs which have been compressed using NewPFD and fibonacci encoding.
 * @param BUF The encoded ECs to decode
 * @param rows The output BUSData array whose EC members are set in this function (up to `row_count`-1).
 * @param row_count The number of elements in `rows` to decode.
 * @param buf_size The size of the input array `BUF`.
 */
void decompress_ec(char *BUF, BUSData *rows, const size_t row_count, const size_t buf_size)
{
	
}

/**
 * @brief Decompress counts which have been compressed using runlen(1) and fibonacci encoding.
 * @param BUF The encoded counts to decode
 * @param rows The output BUSData array whose count members are set in this function (up to `row_count`-1).
 * @param row_count The number of elements in `rows` to decode.
 * @param buf_size The size of the input array `BUF`.
 */
void decompress_counts(char *BUF, BUSData *rows, const size_t row_count, const size_t buf_size) {
	uint64_t *BUF64 = (uint64_t *)BUF;
	size_t buf64_size = (buf_size - 1) / (sizeof(uint64_t)) + 1;

	uint32_t bitpos{0},
		i_fibo{0},
		buf_offset{0};
	size_t row_index = 0;

	uint32_t curr_el = 0,
			 runlen = 0;
	const uint32_t RLE_val{1U};

	while(row_index < row_count){
		curr_el = (uint32_t)fiboDecodeSingle(BUF64, buf64_size, i_fibo, bitpos, buf_offset);
		runlen = curr_el == RLE_val ? (uint32_t)fiboDecodeSingle(BUF64, buf64_size, i_fibo, bitpos, buf_offset) : 1;
		for (int i = 0; i < runlen; ++i){
			rows[row_index].count = curr_el;
			++row_index;
		}
	}
}

/**
 * @brief Decompress counts which have been compressed using runlen(0) and fibonacci encoding.
 * @param BUF The encoded counts to decode
 * @param rows The output BUSData array whose count members are set in this function (up to `row_count`-1).
 * @param row_count The number of elements in `rows` to decode.
 * @param buf_size The size of the input array `BUF`.
 */
void decompress_flags(char *BUF, BUSData *rows, const size_t row_count, const size_t buf_size) {
	uint64_t *BUF64 = (uint64_t *)BUF;
	size_t buf64_size = (buf_size - 1) / (sizeof(uint64_t)) + 1;

	uint32_t bitpos{0},
		i_fibo{0},
		buf_offset{0};
	size_t row_index = 0;

	uint32_t curr_el = 0,
			 runlen = 0;
	const uint32_t RLE_val{0U};

	while(row_index < row_count){
		curr_el = (uint32_t)fiboDecodeSingle(BUF64, buf64_size, i_fibo, bitpos, buf_offset) - 1;
		runlen = curr_el == RLE_val ? (uint32_t)fiboDecodeSingle(BUF64, buf64_size, i_fibo, bitpos, buf_offset) : 1;
		for (int i = 0; i < runlen; ++i)
		{
			rows[row_index].flags = curr_el;
			++row_index;
		}
	}
}


void bustools_decompress(const Bustools_opt &opt)
{
	
}