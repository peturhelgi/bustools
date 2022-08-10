#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <vector>
#include <array>
#include <stdint.h>
#include <zlib.h>

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
template <typename T>
uint64_t fiboDecodeSingle(T const *const buf, const size_t n_buf, uint32_t &i_fibo, uint32_t &bitpos, uint32_t &bufpos)
{

	size_t SIZE = sizeof(T) * 8;
	size_t buf_offset = bufpos;
	size_t bit_offset = SIZE - (bitpos % SIZE) - 1;
	uint64_t num{0};

	int bit = (buf[buf_offset] >> bit_offset) & 1;
	int last_bit = 0;

	++bitpos;

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
void updatePFD(int32_t *primary, std::vector<int32_t> &index_gaps, std::vector<int32_t> &exceptions, uint32_t b_bits, const int32_t min_element)
{
	int i_exception = 0;
	int exception_pos = 0;
	assert(index_gaps.size() == exceptions.size());
	int32_t index = 0;
	int n = index_gaps.size();
	for (int i = 0; i < n; ++i)
	{
		auto index_diff = index_gaps[i];
		uint32_t ex = exceptions[i];
		index += index_diff;
		primary[index] |= (exceptions[i] << b_bits);
	}
	for (int i = 0; i < 512; ++i)
	{
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
size_t PfdParsePrimaryBlock(
	PFD_t *buf,
	const size_t max_elem,
	const int n_ints,
	const uint32_t b_bits,
	int32_t *primary,
	uint32_t &bit_pos,
	size_t buf_offset)
{
	int i = 0;
	constexpr uint32_t wordsize = sizeof(PFD_t) * 8;
	constexpr PFD_t ONE{1};
	while (i < n_ints && buf_offset < max_elem)
	{
		// I: 0 <= bit_pos < 64

		uint32_t n_partial = std::min(b_bits, wordsize - bit_pos);
		uint32_t bits_rem = b_bits - n_partial;

		uint32_t shift = (wordsize - n_partial - bit_pos);
		PFD_t mask = (ONE << n_partial) - 1;

		int32_t num = ((buf[buf_offset] >> shift) & mask) << bits_rem;

		bit_pos = (bit_pos + n_partial) % wordsize;
		// increment buf_offset if bit_pos == 0, i.e. the next element in `buf` is reached.
		buf_offset += (!bit_pos);

		if (bits_rem)
		{
			shift = wordsize - bits_rem;
			num |= buf[buf_offset] >> shift;
			bit_pos += bits_rem;
		}
		primary[i++] = num;
	}

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
void decompress_barcode(char *BUF, BUSData *rows, const size_t row_count, size_t &buf_size)
{
	FIBO_t *fibonacci_buf = (FIBO_t *)BUF;
	size_t fibonacci_bufsize = (buf_size - 1) / (sizeof(FIBO_t)) + 1;

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
		diff = fiboDecodeSingle(fibonacci_buf, fibonacci_bufsize, i_fibo, bitpos, buf_offset) - 1;

		if (diff == RLE_VAL)
		{
			// Runlength decoding
			runlen = fiboDecodeSingle(fibonacci_buf, fibonacci_bufsize, i_fibo, bitpos, buf_offset);
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

	buf_size = (buf_offset + (bitpos > 0)) * sizeof(FIBO_t);
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
void decompress_lossless_umi(char *BUF, BUSData *rows, const size_t row_count, size_t &buf_size) {
	FIBO_t *fibonacci_buf = (FIBO_t *)BUF;
	size_t fibonacci_bufsize = (buf_size - 1) / (sizeof(FIBO_t)) + 1,
		   row_index = 0;

	uint32_t bitpos{0},
		i_fibo{0},
		buf_offset{0};
	uint64_t diff = 0,
			 last_barcode = rows[0].barcode + 1,
			 umi = 0,
			 barcode,
			 runlen = 0;

	const uint64_t RLE_VAL{0ULL};

	while(row_index < row_count){
		diff = fiboDecodeSingle(fibonacci_buf, fibonacci_bufsize, i_fibo, bitpos, buf_offset) - 1;
		barcode = rows[row_index].barcode;
		if (barcode != last_barcode)
		{
			umi = 0;
		}

		if (diff == RLE_VAL)
		{
			// diff is runlen encoded, next values are identical.
			runlen = fiboDecodeSingle(fibonacci_buf, fibonacci_bufsize, i_fibo, bitpos, buf_offset);
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
	buf_size = (buf_offset + (bitpos > 0)) * sizeof(FIBO_t);
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
void decompress_lossy_umi(char *BUF, BUSData *rows, const size_t row_count, size_t &buf_size) {}

/**
 * @brief Decompress ECs which have been compressed using NewPFD and fibonacci encoding.
 * @param BUF The encoded ECs to decode
 * @param rows The output BUSData array whose EC members are set in this function (up to `row_count`-1).
 * @param row_count The number of elements in `rows` to decode.
 * @param buf_size The size of the input array `BUF`.
 */
void decompress_ec(char *BUF, BUSData *rows, const size_t row_count, size_t &buf_size)
{
	uint32_t b_bits = 1;
	uint32_t buf_offset{0};

	PFD_t *pfd_buf = (PFD_t *)BUF;

	size_t pfd_bufsize = (buf_size - 1) / sizeof(PFD_t) + 1;
	uint32_t bitpos{0}, i_fibo{0};

	int32_t min_element{0};
	constexpr size_t block_size{512};

	int32_t primary[block_size];
	size_t n_exceptions{0};
	std::vector<int32_t> exceptions;
	std::vector<int32_t> index_gaps;

	b_bits = fiboDecodeSingle(pfd_buf, pfd_bufsize, i_fibo, bitpos, buf_offset) - 1;
	min_element = (int32_t)fiboDecodeSingle(pfd_buf, pfd_bufsize, i_fibo, bitpos, buf_offset) - 1;
	n_exceptions = fiboDecodeSingle(pfd_buf, pfd_bufsize, i_fibo, bitpos, buf_offset) - 1;

	size_t row_index = 0;
	while (b_bits)
	{
		index_gaps.clear();
		exceptions.clear();

		for (int i = 0; i < n_exceptions; ++i)
		{
			index_gaps.push_back(fiboDecodeSingle(pfd_buf, pfd_bufsize, i_fibo, bitpos, buf_offset) - 1);
		}
		
		for (int i = 0; i < n_exceptions; ++i)
		{
			exceptions.push_back(fiboDecodeSingle(pfd_buf, pfd_bufsize, i_fibo, bitpos, buf_offset));
		}

		buf_offset += (bitpos > 0);
		bitpos = 0;

		buf_offset = PfdParsePrimaryBlock(pfd_buf, pfd_bufsize, 512, b_bits, primary, bitpos, buf_offset);
		updatePFD(primary, index_gaps, exceptions, b_bits, min_element);
		for (int i = 0; i < block_size && row_index < row_count; ++i){
			rows[row_index].ec = primary[i];
			++row_index;
		}
		
		b_bits = fiboDecodeSingle(pfd_buf, pfd_bufsize, i_fibo, bitpos, buf_offset) - 1;
		min_element = (int32_t)fiboDecodeSingle(pfd_buf, pfd_bufsize, i_fibo, bitpos, buf_offset) - 1;
		n_exceptions = fiboDecodeSingle(pfd_buf, pfd_bufsize, i_fibo, bitpos, buf_offset) - 1;
	}
	buf_size = (buf_offset + (bitpos > 0)) * sizeof(PFD_t);
}

/**
 * @brief Decompress counts which have been compressed using runlen(1) and fibonacci encoding.
 * @param BUF The encoded counts to decode
 * @param rows The output BUSData array whose count members are set in this function (up to `row_count`-1).
 * @param row_count The number of elements in `rows` to decode.
 * @param buf_size The size of the input array `BUF`.
 */
void decompress_counts(char *BUF, BUSData *rows, const size_t row_count, size_t &buf_size) {
	FIBO_t *fibonacci_buf = (FIBO_t *)BUF;
	size_t fibonacci_bufsize = (buf_size - 1) / (sizeof(FIBO_t)) + 1;

	uint32_t bitpos{0},
		i_fibo{0},
		buf_offset{0};
	size_t row_index = 0;

	uint32_t curr_el = 0,
			 runlen = 0;
	const uint32_t RLE_val{1U};

	while(row_index < row_count)
	{
		curr_el = (uint32_t)fiboDecodeSingle(fibonacci_buf, fibonacci_bufsize, i_fibo, bitpos, buf_offset);
		runlen = curr_el == RLE_val ? (uint32_t)fiboDecodeSingle(fibonacci_buf, fibonacci_bufsize, i_fibo, bitpos, buf_offset) : 1;
		for (int i = 0; i < runlen; ++i){
			rows[row_index].count = curr_el;
			++row_index;
		}
	}
	buf_size =(buf_offset + (bitpos > 0)) * sizeof(FIBO_t);
}

/**
 * @brief Decompress counts which have been compressed using runlen(0) and fibonacci encoding.
 * @param BUF The encoded counts to decode
 * @param rows The output BUSData array whose count members are set in this function (up to `row_count`-1).
 * @param row_count The number of elements in `rows` to decode.
 * @param buf_size The size of the input array `BUF`.
 */
void decompress_flags(char *BUF, BUSData *rows, const size_t row_count, size_t &buf_size) {
	FIBO_t *fibonacci_buf = (FIBO_t *)BUF;
	size_t fibonacci_bufsize = (buf_size - 1) / (sizeof(FIBO_t)) + 1,
		row_index = 0;

	uint32_t bitpos{0},
		i_fibo{0},
		buf_offset{0},
		curr_el{0},
		runlen{0};

	const uint32_t RLE_val{0U};

	while(row_index < row_count)
	{
		curr_el = (uint32_t)fiboDecodeSingle(fibonacci_buf, fibonacci_bufsize, i_fibo, bitpos, buf_offset) - 1;
		runlen = curr_el == RLE_val ? (uint32_t)fiboDecodeSingle(fibonacci_buf, fibonacci_bufsize, i_fibo, bitpos, buf_offset) : 1;
		for (int i = 0; i < runlen; ++i)
		{
			rows[row_index].flags = curr_el;
			++row_index;
		}
	}
	buf_size = (buf_offset + (bitpos > 0)) * sizeof(FIBO_t);
	
}

template <typename T>
void decompress_barcode_fibo(char *BUF, BUSData *rows, const size_t row_count, size_t &buf_size)
{
	T *fibonacci_buf = (T *)BUF;
	size_t fibonacci_bufsize = (buf_size - 1) / sizeof(T) + 1;
	uint32_t bitpos{0},
		i_fibo{0},
		buf_offset{0};
	
	for (int i = 0; i < row_count; ++i){
		rows[i].barcode = fiboDecodeSingle(fibonacci_buf, fibonacci_bufsize, i_fibo, bitpos, buf_offset) - 1;
	}
	buf_size = (buf_offset + (bitpos > 0)) * sizeof(T);
}
template <typename T>
void decompress_UMI_fibo(char *BUF, BUSData *rows, const size_t row_count, size_t &buf_size)
{
	T *fibonacci_buf = (T *)BUF;
	size_t fibonacci_bufsize = (buf_size - 1) / sizeof(T) + 1;
	uint32_t bitpos{0},
		i_fibo{0},
		buf_offset{0};

	for (int i = 0; i < row_count; ++i)
	{
		rows[i].UMI = fiboDecodeSingle(fibonacci_buf, fibonacci_bufsize, i_fibo, bitpos, buf_offset) - 1;
	}

	buf_size = (buf_offset + (bitpos > 0)) * sizeof(T);
}
template <typename T>
void decompress_EC_fibo(char *BUF, BUSData *rows, const size_t row_count, size_t &buf_size)
{
	T *fibonacci_buf = (T *)BUF;
	size_t fibonacci_bufsize = (buf_size - 1) / sizeof(T) + 1;
	uint32_t bitpos{0},
		i_fibo{0},
		buf_offset{0};

	for (int i = 0; i < row_count; ++i)
	{
		rows[i].ec = fiboDecodeSingle(fibonacci_buf, fibonacci_bufsize, i_fibo, bitpos, buf_offset) - 1;
	}

	buf_size = (buf_offset + (bitpos > 0)) * sizeof(T);
}
template <typename T>
void decompress_count_fibo(char *BUF, BUSData *rows, const size_t row_count, size_t &buf_size)
{
	T *fibonacci_buf = (T *)BUF;
	size_t fibonacci_bufsize = (buf_size - 1) / sizeof(T) + 1;
	uint32_t bitpos{0},
		i_fibo{0},
		buf_offset{0};

	for (int i = 0; i < row_count; ++i)
	{
		rows[i].count = fiboDecodeSingle(fibonacci_buf, fibonacci_bufsize, i_fibo, bitpos, buf_offset);
	}

	buf_size = (buf_offset + (bitpos > 0)) * sizeof(T);
}
template <typename T>
void decompress_flags_fibo(char *BUF, BUSData *rows, const size_t row_count, size_t &buf_size)
{
	T *fibonacci_buf = (T *)BUF;
	size_t fibonacci_bufsize = (buf_size - 1) / sizeof(T) + 1;
	uint32_t bitpos{0},
		i_fibo{0},
		buf_offset{0};

	for (int i = 0; i < row_count; ++i)
	{
		rows[i].flags = fiboDecodeSingle(fibonacci_buf, fibonacci_bufsize, i_fibo, bitpos, buf_offset) - 1;
	}

	buf_size = (buf_offset + (bitpos > 0)) * sizeof(T);
}

void decompress_barcode_zlib(char *BUF, BUSData *rows, const size_t row_count, size_t &buf_size)
{
	uint64_t *buf = new uint64_t[row_count];
	uLongf dest_len = row_count * sizeof(uint64_t);
	uLong src_len = buf_size;
	int status = uncompress2((Bytef *)buf, &dest_len, (Bytef *)BUF, &src_len);
	if (status != Z_OK)
	{
		std::cerr << "zlib error: " << status << "\n";
	}
	for (int i = 0; i < row_count; ++i)
	{
		rows[i].barcode = buf[i];
	}
	buf_size = src_len;
}


void decompress_UMI_zlib(char *BUF, BUSData *rows, const size_t row_count, size_t &buf_size)
{
	uint64_t *buf = new uint64_t[row_count];
	uLongf dest_len = row_count * sizeof(uint64_t);
	uLong src_len = buf_size;
	int status = uncompress2((Bytef *)buf, &dest_len, (Bytef *)BUF, &src_len);
	if (status != Z_OK)
	{
		std::cerr << "zlib error: " << status << "\n";
	}
	for (int i = 0; i < row_count; ++i)
	{
		rows[i].UMI = buf[i];
	}
	buf_size = src_len;
}

void decompress_EC_zlib(char *BUF, BUSData *rows, const size_t row_count, size_t &buf_size)
{
	int32_t *buf = new int32_t[row_count];
	uLongf dest_len = row_count * sizeof(int32_t);
	uLong src_len = buf_size;
	int status = uncompress2((Bytef *)buf, &dest_len, (Bytef *)BUF, &src_len);
	if (status != Z_OK)
	{
		std::cerr << "zlib error: " << status << "\n";
	}
	for (int i = 0; i < row_count; ++i)
	{
		rows[i].ec = buf[i];
	}
	buf_size = src_len;
}

void decompress_count_zlib(char *BUF, BUSData *rows, const size_t row_count, size_t &buf_size)
{
	uint32_t *buf = new uint32_t[row_count];
	uLongf dest_len = row_count * sizeof(uint32_t);
	uLong src_len = buf_size;
	int status = uncompress2((Bytef *)buf, &dest_len, (Bytef *)BUF, &src_len);
	if (status != Z_OK)
	{
		std::cerr << "zlib error: " << status << "\n";
	}
	for (int i = 0; i < row_count; ++i)
	{
		rows[i].count = buf[i];
	}
	buf_size = src_len;
}

void decompress_flags_zlib(char *BUF, BUSData *rows, const size_t row_count, size_t &buf_size)
{
	uint32_t *buf = new uint32_t[row_count];
	uLongf dest_len = row_count * sizeof(uint32_t);
	uLong src_len = buf_size;
	int status = uncompress2((Bytef *)buf, &dest_len, (Bytef *)BUF, &src_len);
	if (status != Z_OK)
	{
		std::cerr << "zlib error: " << status << "\n";
	}
	for (int i = 0; i < row_count; ++i)
	{
		rows[i].flags = buf[i];
	}
	buf_size = src_len;
}

typedef void (*decompress_ptr)(char *, BUSData *, const size_t, size_t &);

void select_decompressors(const compressed_BUSHeader &comp_h, decompress_ptr decompressors[5]){
	decompress_ptr fibo[5]{
		&decompress_barcode_fibo<FIBO_t>,
		&decompress_UMI_fibo<FIBO_t>,
		&decompress_EC_fibo<FIBO_t>,
		&decompress_count_fibo<FIBO_t>,
		&decompress_flags_fibo<FIBO_t>,
	};
	decompress_ptr zlib[5]{
		&decompress_barcode_zlib,
		&decompress_UMI_zlib,
		&decompress_EC_zlib,
		&decompress_count_zlib,
		&decompress_flags_zlib,
	};

	for (int i = 0; i < 5; ++i)
	{
		if(comp_h.fibo_zlib_compress & (1 << (5 + 4-i))){
			decompressors[i] = fibo[i];
		}
		else if (comp_h.fibo_zlib_compress & (1 << (4 - i))){
			decompressors[i] = zlib[i];
		}
	}
}

void decompress_block(
	const size_t row_count,
	const size_t bufsize,
	const decompress_ptr decompressors[5],
	char *const BUF,
	BUSData *busdata)
{
	size_t buf_pos = 0;
	for (int i = 0; i < 5; ++i)
	{
		size_t srclen = bufsize;
		decompressors[i](BUF + buf_pos, busdata, row_count, srclen);
		buf_pos += srclen;
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
		std::cerr << "Error: Failed to parse header.\n";
		return;
	}

	std::vector<uint32_t> chunk_sizes(comp_header.n_chunks+1);
	in.read((char*)&chunk_sizes[0], (comp_header.n_chunks + 1) * sizeof(uint32_t));

	writeHeader(outf, comp_header.extra_header);
	BUSData *busdata = new BUSData[comp_header.chunk_size];

	decompress_ptr decompress_umi = comp_header.lossy_umi ? &decompress_lossy_umi : &decompress_lossless_umi;
	decompress_ptr decompressors[5]{
		&decompress_barcode,
		decompress_umi,
		&decompress_ec,
		&decompress_counts,
		&decompress_flags,
	};
	select_decompressors(comp_header, decompressors);

	uint32_t curr_row_count = comp_header.chunk_size,
			 i_chunk = 0;

	const auto chunk_size_end = chunk_sizes.end();
	try
	{
		uint32_t max_block_size = *std::max_element(chunk_sizes.begin(), chunk_sizes.end());
		char *BUF = new char[max_block_size];
		for (auto chunk_size_it = chunk_sizes.begin(); in.good() && chunk_size_it < chunk_size_end; ++chunk_size_it, ++i_chunk)
		{
			in.read((char *)BUF, *chunk_size_it);

			curr_row_count = (i_chunk == comp_header.n_chunks) ? comp_header.last_chunk : curr_row_count;
			decompress_block(curr_row_count, *chunk_size_it, decompressors, BUF, busdata);
			outf.write((char *)busdata, curr_row_count * sizeof(BUSData));
		}

		delete[] BUF;
	}
	catch (std::bad_alloc)
	{
		std::cerr << "Unable to allocate bytes" << std::endl;
	}

	delete[] busdata;
}