#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <vector>
#include <array>
#include <stdint.h>

#include "Common.hpp"
#include "BUSData.h"
#include "bustools_compress.h"

/**
 * @brief Encode `num` using fibonacci encoding into buf, starting at bitpos.
 * @pre the next 128 bits in `buf` are 0 starting from `bitpos`, and wrapped around 192.
 * 		0 <= bit_pos < 3*64 == 192
 *
 * @pre num > 0
 * @post buf now contains the fibo encoding of num from [pre(bitpos); post(bitpos)].
 *		bit_pos has changed, 0 <= bitpos < 3*64 == 192.
 * 		0 or more elements in buf have been concentrated with fibo-encoded numbers and thus written to obuf.
 *
 * @note The largest fibonacci number that fits in a uint64_t is the 91st (0-indexed) fibo number, i.e. the first 92 fibonacci numbers fit in a 64-bit uint.
 * 	Fibonacci encoding appends a single bit as a stop code. Hence, the longest fibonacci encoding uses 93 bits.
 *
 * @param num the number to encode, num > 0
 * @param buf array of 3 uint64_t. The encoded bits are stored in buf.
 * @param bitpos the bit position in buf where the fibo encoding of num starts.
 * @param obuf the ostream buffer to write to.
 */
void fiboEncode(const uint64_t num, uint64_t buf[3], uint32_t &bitpos, std::ostream &obuf)
{
	constexpr uint32_t max_bitpos = 3 * 64;
	const uint32_t curr_byte_pos = bitpos / 64;

	const auto fibo_begin = fibo64.begin();
	const auto fibo_end = fibo64.end();

	uint64_t remainder = num;

	// the ith fibonacci number is the largest fibo not greater than remainder
	auto i = std::upper_bound(fibo_begin, fibo_end, remainder) - 1;

	const uint32_t n_bits = (i - fibo_begin) + 2;
	uint32_t next_bit_pos = bitpos + n_bits - 1,
			 bit_offset = next_bit_pos % 64,
			 buf_offset = (next_bit_pos / 64) % 3;

	// Set the stop bit.
	buf[buf_offset] |= 1ULL << (63 - bit_offset);

	i = fibo_end;
	while (remainder > 0)
	{
		i = std::upper_bound(fibo_begin, i, remainder) - 1;
		next_bit_pos = bitpos + (i - fibo_begin);
		buf_offset = (next_bit_pos / 64) % 3;
		bit_offset = next_bit_pos % 64;

		buf[buf_offset] |= 1ULL << (63 - bit_offset);
		remainder -= *i;
	}

	// n_elems is the number of saturated elements in buf.
	int n_elems = (bitpos + n_bits) / 64 - curr_byte_pos;
	bitpos = (bitpos + n_bits) % max_bitpos;

	// write fibo-encoded values to output buffer for concentrated elements in buf.
	for (int p = 0; p < n_elems; ++p)
	{
		auto &elem = buf[(curr_byte_pos + p) % 3];
		obuf.write((char *)&elem, sizeof(elem));
		elem = 0;
	}
}

/**
 * @brief Compress barcodes of rows using delta-runlen(0)-fibonacci encoding and write to `of`.
 *
 * @param rows BUSData array, contains at least `row_count` elements
 * @param row_count The number of barcodes to compress.
 * @param of The ostream for writing the encoding to.
 */
void compress_barcodes(BUSData const *const rows, const int row_count, std::ostream &of)
{
	uint64_t barcode, last_bc = 0;
	uint64_t runlen = 0;

	uint64_t buf[3]{0, 0, 0};
	uint32_t bit_pos{0};

	for (int i = 0; i < row_count; ++i)
	{
		barcode = rows[i].barcode;

		// delta encoding
		barcode -= last_bc;

		// Runlength encoding of zeros
		if (barcode == 0) {
			++runlen;
		} else
		{
			// Increment values as fibo cannot encode 0
			if (runlen) {
				fiboEncode(1ULL, buf, bit_pos, of);
				fiboEncode(runlen, buf, bit_pos, of);
				runlen = 0;
			}
			fiboEncode(barcode + 1, buf, bit_pos, of);
		}
		last_bc = rows[i].barcode;
	}

	// Take care of the last run of zeros in the delta-encoded barcodes
	if (runlen) {
		fiboEncode(1ULL, buf, bit_pos, of);
		fiboEncode(runlen, buf, bit_pos, of);
	}

	// Write out the last fibonacci element if applicable.
	if (bit_pos % 64) {
		auto &element = buf[bit_pos / 64];
		of.write((char *)(&element), sizeof(element));
	}
}

void lossless_compress_umis(BUSData const *const rows, const int row_count, std::ostream &of)
{
	uint64_t last_bc = 0,
			 last_umi = 0,
			 bc, umi, diff;

	uint64_t buf[3]{0, 0, 0};
	uint32_t bit_pos{0};
	const uint32_t RLE_val{0ULL};
	uint64_t runlen{0};

	for (int i = 0; i < row_count; ++i){
		bc = rows[i].barcode;
		umi = rows[i].UMI;

		if(last_bc != bc){
			last_umi = 0;
		}
		diff = umi - last_umi;
		last_umi = umi;
		last_bc = bc;

		if(diff == RLE_val){
			++runlen;
		}
		else{
			if(runlen){
				fiboEncode(RLE_val + 1, buf, bit_pos, of);
				fiboEncode(runlen, buf, bit_pos, of);
				runlen = 0;
			}
			fiboEncode(diff + 1, buf, bit_pos, of);
		}
	}
	if(runlen){
		fiboEncode(RLE_val + 1, buf, bit_pos, of);
		fiboEncode(runlen, buf, bit_pos, of);
	}

	// Write last bytes when the buffers have not been completely filled.
	if(bit_pos % 64)
	{
		auto &element = buf[bit_pos / 64];
		of.write((char *)(&element), sizeof(element));
	}
}

void lossy_compress_umis(BUSData const *const rows, const int row_count, std::ostream &of)
{
	std::cout << "BEWARE: Lossy compression\n";

}

void compress_ecs(BUSData const *const rows, const int row_count, std::ostream &of)
{
}

/**
 * @brief Compress counts of rows using runlength(1)-fibonacci encoding and write to `of`.
 *
 * @param rows BUSData array, contains at least `row_count` elements
 * @param row_count The number of counts to compress.
 * @param of The ostream for writing the encoding to.
 */
void compress_counts(BUSData const *const rows, const int row_count, std::ostream &of)
{
	const uint32_t RLE_val{1UL};
	uint32_t count,
		runlen{0},
		bitpos{0};

	uint64_t buf[3]{0, 0, 0};

	for (int i = 0; i < row_count; ++i)
	{
		count = rows[i].count;

		if (count == RLE_val)
		{
			++runlen;
		}
		else
		{
			if (runlen)
			{
				// Runlength-encode 1s.
				fiboEncode(RLE_val, buf, bitpos, of);
				fiboEncode(runlen, buf, bitpos, of);
				runlen = 0;
			}

			fiboEncode(count, buf, bitpos, of);
		}
	}
	if (runlen)
	{
		// Runlength-encode last run of 1s.
		fiboEncode(RLE_val, buf, bitpos, of);
		fiboEncode(runlen, buf, bitpos, of);
	}

	// Write last bytes when if the fibonacci_buf has not been saturated.
	if (bitpos % 64)
	{
		auto &element = buf[bitpos / 64];
		of.write((char *)(&element), sizeof(element));
	}
}

void compress_flags(BUSData const *const rows, const int row_count, std::ostream &of)
{
	uint32_t flag,
		runlen = 0;

	uint64_t buf[3]{0, 0, 0};
	uint32_t bit_pos{0};

	const uint32_t RLE_val{0UL};

	for (int i = 0; i < row_count; ++i)
	{
		flag = rows[i].flags;

		// Runlength, encode 1
		if (flag == RLE_val)
		{
			++runlen;
		}
		else
		{
			// Increment values as fibo cannot encode 0
			if (runlen)
			{
				// fibonacci encode 0(+1) and runlen
				fiboEncode(RLE_val + 1, buf, bit_pos, of);
				fiboEncode(runlen, buf, bit_pos, of);
				runlen = 0;
			}

			// encode value + 1
			fiboEncode(flag + 1, buf, bit_pos, of);
		}
	}
	if (runlen)
	{
		fiboEncode(RLE_val + 1, buf, bit_pos, of);
		fiboEncode(runlen, buf, bit_pos, of);
	}

	// TODO: check if we need to change this.
	if (bit_pos % 64)
	{
		auto &element = buf[bit_pos / 64];
		of.write((char *)(&element), sizeof(element));
	}
}

void bustools_compress(const Bustools_opt &opt){
	BUSHeader h;
	
	// maximum chunk_size set by max_memory
	const size_t ROW_SIZE = sizeof(BUSData);
	size_t N = opt.max_memory / ROW_SIZE;
	const size_t chunk_size = (N < opt.chunk_size) ? N : opt.chunk_size;

	void (*funcs [])(BUSData const *const, const int, std::ostream &) = {&lossless_compress_umis, &lossy_compress_umis};
	if(opt.lossy_umi){
		std::cerr << "Lossy UMI has not been implemented, using lossless instead" << std::endl;
	}
	void (*compress_umis)(BUSData const *const, const int, std::ostream&) = funcs[0];

	BUSData *p = new BUSData[chunk_size];

	std::ofstream of;
	std::streambuf *buf = nullptr;

	if (opt.stream_out){
		buf = std::cout.rdbuf();
	}
	else{
		of.open(opt.output);
		buf = of.rdbuf();
	}

	std::ostream outf(buf);

	for (const auto &infn : opt.files)
	{
		std::streambuf *inbuf;
		std::ifstream inf;
		if(opt.stream_in){
			inbuf = std::cin.rdbuf();
		}
		else{
			inf.open(infn.c_str(), std::ios::binary);
			inbuf = inf.rdbuf();
		}
		std::istream in(inbuf);

		parseHeader(in, h);
		compressed_BUSHeader comp_h;
		comp_h.chunk_size = chunk_size;
		comp_h.lossy_umi = false; //opt.lossy_umi;
		// copy header to comp_h.extra_header;

		uint64_t block_counter = 0;
		size_t last_row_count = 0;
		while (in.good())
		{
			in.read((char *)p, chunk_size * ROW_SIZE);
			size_t row_count = in.gcount() / ROW_SIZE;
			last_row_count = row_count;
			++block_counter;

			compress_barcodes(p, row_count, outf);
			compress_umis(p, row_count, outf);
			// compress_ecs(p, row_count, outf);
			compress_counts(p, row_count, outf);
			compress_flags(p, row_count, outf);
		}

		--block_counter;
		comp_h.last_chunk = last_row_count;
		comp_h.n_chunks = block_counter;
	}
}