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
 * @brief Finish NewPFD encoding after parsing to account for exceptions.
 * 
 * @param primary A pointer to a block of packed fixed-width integers each of size `b_bits`.
 * @param index_gaps Delta encoded indices of exceptions in `primary`.
 * @param exceptions The most significant bits of exceptions, each one is always > 0.
 * @param b_bits The number of bits used for packing numbers into `primary`.
 */
void updatePFD(int32_t *primary, std::vector<int32_t> &index_gaps, std::vector<int32_t> &exceptions, uint32_t b_bits, const int32_t min_element){
	int i_exception = 0;
	int exception_pos = 0;
	assert(index_gaps.size() == exceptions.size());
	int32_t index = 0;
	int n = index_gaps.size();
	for (int i = 0; i < n; ++i){
		auto index_diff = index_gaps[i];
		uint32_t ex = exceptions[i];
		index += index_diff;
		primary[index] |= (exceptions[i] << b_bits);
	}
	for (int i = 0; i < 512; ++i){
		primary[i] += min_element;
	}
}

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
	uint32_t b_bits = 1;
	uint32_t buf_offset{0};

	uint64_t *BUF64 = (uint64_t *)BUF;

	size_t buf64_size = (buf_size - 1) / sizeof(uint64_t) + 1;
	uint32_t bitpos{0}, i_fibo{0};

	int32_t min_element{0};

	int32_t primary[512];
	size_t n_exceptions{0};
	std::vector<int32_t> exceptions;
	std::vector<int32_t> index_gaps;

	b_bits = fiboDecodeSingle(BUF64, buf64_size, i_fibo, bitpos, buf_offset) - 1;
	min_element = (int32_t)fiboDecodeSingle(BUF64, buf64_size, i_fibo, bitpos, buf_offset) - 1;
	n_exceptions = fiboDecodeSingle(BUF64, buf64_size, i_fibo, bitpos, buf_offset) - 1;

	size_t row_index = 0;
	while (b_bits)
	{

		index_gaps.clear();
		exceptions.clear();

		const auto *start_pos = BUF;

		for (int i = 0; i < n_exceptions; ++i)
		{
			index_gaps.push_back(fiboDecodeSingle(BUF64, buf64_size, i_fibo, bitpos, buf_offset) - 1);
		}
		
		for (int i = 0; i < n_exceptions; ++i)
		{
			exceptions.push_back(fiboDecodeSingle(BUF64, buf64_size, i_fibo, bitpos, buf_offset));
		}

		buf_offset += (bitpos > 0);
		bitpos = 0;

		buf_offset = PfdParsePrimaryBlock(BUF64, buf64_size, 512, b_bits, primary, bitpos, buf_offset);
		updatePFD(primary, index_gaps, exceptions, b_bits, min_element);
		for (int i = 0; i < 512 && row_index < row_count; ++i){
			rows[row_index].ec = primary[i];
			++row_index;
		}
		
		b_bits = fiboDecodeSingle(BUF64, buf64_size, i_fibo, bitpos, buf_offset) - 1;
		min_element = (int32_t)fiboDecodeSingle(BUF64, buf64_size, i_fibo, bitpos, buf_offset) - 1;
		n_exceptions = fiboDecodeSingle(BUF64, buf64_size, i_fibo, bitpos, buf_offset) - 1;
	}
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
	compressed_BUSHeader comp_header;
	std::ofstream of;
	std::streambuf *buf = nullptr;

	if (opt.stream_out)
	{
		buf = std::cout.rdbuf();
	}
	else
	{
		of.open(opt.output);
		buf = of.rdbuf();
	}

	std::ostream outf(buf);
	const auto &infn = opt.files[0];
	std::streambuf *inbuf;
	std::ifstream inf;

	// don't allow stream in for now.
	const bool allow_streamin = false;
	if (opt.stream_in && allow_streamin)
	{
		inbuf = std::cin.rdbuf();
	}
	else
	{
		inf.open(infn.c_str(), std::ios::binary);
		inbuf = inf.rdbuf();
	}

	std::istream in(inbuf);

	if(!parseCompressedHeader(in, comp_header)){
		std::cerr << "Error: Failed at parsing header.\n";
		return;
	}

	writeHeader(outf, comp_header.extra_header);
	BUSData *busdata = new BUSData[comp_header.chunk_size];


	const auto data_pos = in.tellg();
	int n_decompressors = 5;
	uint32_t n_cols = n_decompressors * (comp_header.n_chunks + 1);
	int32_t offset = n_cols * sizeof(uint32_t);

	in.seekg(-offset, std::ios_base::end);

	std::vector<uint32_t> col_sizes(n_cols);
	in.read((char *)&col_sizes[0], offset);


	auto col_size_it = col_sizes.begin();
	const auto col_size_end = col_sizes.end();

	void (*decompress_umi)(char *, BUSData *, const size_t, const size_t) = comp_header.lossy_umi ? &decompress_lossy_umi : &decompress_lossless_umi;

	void (*decompressors[])(char *, BUSData *, const size_t, const size_t) = {
		&decompress_barcode,
		decompress_umi,
		&decompress_ec,
		&decompress_counts,
		&decompress_flags,
	};

	

	in.seekg(data_pos, std::ios_base::beg);

	uint32_t curr_chunk_size = comp_header.chunk_size;
	uint64_t i_chunk = 0;
	size_t total = 0;
	int i_decompressor = 0;

	while (in.good() && col_size_it < col_size_end)
	{
		try{
			char *BUF = new char[*col_size_it];
			curr_chunk_size = (i_chunk == comp_header.n_chunks) ? comp_header.last_chunk : curr_chunk_size;

			in.read((char *)&BUF[0], *col_size_it);

			decompressors[i_decompressor](BUF, busdata, curr_chunk_size, *col_size_it);

			i_decompressor = (i_decompressor + 1) % n_decompressors;
			i_chunk += !i_decompressor;
			++col_size_it;

			if (i_decompressor == 0)
				outf.write((char *)busdata, curr_chunk_size * sizeof(BUSData));

			delete[] BUF;
		}
		catch(std::bad_alloc){
			std::cerr << "Unable to allocate bytes" << std::endl;
		}
	}

	delete[] busdata;
}