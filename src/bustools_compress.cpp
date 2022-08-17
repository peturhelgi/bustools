#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <vector>
#include <array>
#include <stdint.h>
#include <zlib.h>

#include "Common.hpp"
#include "BUSData.h"
#include "bustools_compress.h"

constexpr size_t cache_size = 256;

template <typename T>
inline void flush_fibonacci(T buf[], const uint32_t &bitpos, std::ostream &of)
{
	uint32_t sz = sizeof(T) * 8;
	if(bitpos % sz){
		auto &element = buf[bitpos / sz];
		of.write((char *)&element, sizeof(element));
	}
}

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
template <size_t bufsize, typename BUF_t>
void fiboEncode(const uint64_t num, BUF_t *buf, uint32_t &bitpos, std::ostream &obuf)
{
	constexpr uint32_t word_size = sizeof(BUF_t) * 8;
	constexpr uint32_t max_bitpos = bufsize * word_size;
	constexpr BUF_t ONE{1};

	const uint32_t curr_byte_pos = bitpos / word_size;

	const auto fibo_begin = fibo64.begin();
	const auto fibo_end = fibo64.end();

	uint64_t remainder = num;

	// the ith fibonacci number is the largest fibo not greater than remainder
	auto i = std::upper_bound(fibo_begin, fibo_end, remainder) - 1;

	const uint32_t n_bits = (i - fibo_begin) + 2;
	uint32_t next_bit_pos = bitpos + n_bits - 1,
			 bit_offset = next_bit_pos % word_size,
			 buf_offset = (next_bit_pos / word_size) % bufsize;

	// Set the stop bit.
	buf[buf_offset] |= shifted_64[word_size - 1 - bit_offset]; // ONE << (word_size - 1 - bit_offset);

	++i;
	while (remainder > 0)
	{
		i = std::upper_bound(fibo_begin, i, remainder) - 1;
		next_bit_pos = bitpos + (i - fibo_begin);
		buf_offset = (next_bit_pos / word_size) % bufsize;
		bit_offset = next_bit_pos % word_size;

		buf[buf_offset] |= shifted_64[word_size - 1 - bit_offset]; // ONE << (word_size - 1 - bit_offset);
		remainder -= *i;
	}

	// n_elems is the number of saturated elements in buf.
	int n_elems = (bitpos + n_bits) / word_size - curr_byte_pos;
	bitpos = (bitpos + n_bits) % max_bitpos;

	// write fibo-encoded values to output buffer for concentrated elements in buf.
	for (int p = 0; p < n_elems; ++p)
	{
		auto &elem = buf[(curr_byte_pos + p) % bufsize];
		obuf.write((char *)&elem, sizeof(elem));
		elem = 0;
	}
}

template <typename T>
uint64_t eliasDeltaDecode(T *&buf, uint32_t &bitpos)
{
	size_t wordsize = sizeof(T) * 8;
	size_t max_pos = wordsize - 1;
	uint32_t L = 0, N = 0;
	uint32_t bit{0};

	// Decode unary encoding of L
	while (!(buf[bitpos / wordsize] >> (max_pos - bitpos % wordsize) & 1))
	{
		++L;
		++bitpos;
	}

	// Decode binary encoding of N+1
	for (int i = 0; i < L + 1; ++i)
	{
		bit = buf[bitpos / wordsize] >> (max_pos - bitpos % wordsize) & 1;
		N += bit << (L - i);
		++bitpos;
	}
	--N;

	uint64_t num = 1ULL << N;
	// Decode headless binary encoding of num
	for (int i = 0; i < N; ++i)
	{
		bit = buf[bitpos / wordsize] >> (max_pos - bitpos % wordsize) & 1;
		num += bit << (N - i - 1);
		++bitpos;
	}

	buf += bitpos / wordsize;
	bitpos %= wordsize;
	return num;
}

template <typename T>
void printout(const T num)
{
	size_t bits = sizeof(num) * 8;
	for (int i = bits - 1; i >= 0; --i)
		std::cout << (num >> i & 1);
	std::cout << '\n';
}

template <typename SRC_T, typename DEST_T>
void pack_num(
	const uint32_t b_bits,
	SRC_T elem,
	DEST_T *buf,
	uint32_t &bitpos)
{
	constexpr size_t wordsize = sizeof(DEST_T) * 8;
	int nbits_left = b_bits;
	int chunk = wordsize - (bitpos % wordsize);
	DEST_T ONE{1};
	DEST_T num{elem};
	nbits_left -= chunk;
	while (nbits_left > 0)
	{
		DEST_T mask = (ONE << chunk) - 1;
		buf[bitpos / wordsize] |= (num >> nbits_left) & mask;
		bitpos += chunk;

		chunk = wordsize - (bitpos % wordsize);
		nbits_left -= chunk;
	}

	buf[bitpos / wordsize] |= num << -nbits_left;
	bitpos += chunk + nbits_left;
}

template <typename BUF_T>
void eliasDeltaEncode(const int i, const uint64_t num, BUF_T *&buf, uint32_t &bitpos)
{
	// Largest num is 2^64 - 1 with largest N+1 == 64
	// Thus, an array L_vals[64] = {31 - clz(n+1) for n in range(64)} precompute Ls;
	constexpr size_t wordsize = sizeof(BUF_T) * 8;

	uint32_t N = 63 - __builtin_clzll(num);
	uint32_t L = 31 - __builtin_clz(N + 1);

	// Unary encoding of L is implicitly added via shifts in the output buffer.
	// Thus, the encoding is binary(N+1) prepended to the headlessBinary(num).
	// We can use a precomputed array for (N + 1) << N, and (1 << N) for 0 <= N < 64.

	bitpos += L; // unary encode L by assuming zeros in buffer
	pack_num(L + 1, N + 1, buf, bitpos);
	pack_num(N, num - (1 << N), buf, bitpos);

	buf += bitpos / wordsize;
	bitpos %= wordsize;
}

template <size_t bufsize>
void fiboEncode8_zero(const uint64_t num, unsigned char buf[24], uint32_t &bitpos, std::ostream &obuf)
{
	constexpr uint32_t word_size = sizeof(unsigned char) * 8;
	constexpr uint32_t max_bitpos = bufsize * word_size;

	const uint32_t curr_byte_pos = bitpos / word_size;

	const auto fibo_begin = fibo64_0.begin();
	const auto fibo_end = fibo64_0.end();

	uint64_t remainder = num;

	// the ith fibonacci number is the largest fibo not greater than remainder
	auto i = std::upper_bound(fibo_begin, fibo_end, remainder) - 1;

	const uint32_t n_bits = (i - fibo_begin) + 2;
	uint32_t next_bit_pos = bitpos + n_bits - 1,
			 bit_offset = next_bit_pos % word_size,
			 buf_offset = (next_bit_pos / word_size) % bufsize;

	// Set the stop bit.
	buf[buf_offset] |= shifted_8[word_size - 1 - bit_offset];

	++i;
	do
	{
		i = std::upper_bound(fibo_begin, i, remainder) - 1;
		next_bit_pos = bitpos + (i - fibo_begin);
		buf_offset = (next_bit_pos / word_size) % bufsize;
		bit_offset = next_bit_pos % word_size;

		buf[buf_offset] |= shifted_8[word_size - 1 - bit_offset];
		remainder -= *i;
	} while (remainder > 0);

	// n_elems is the number of saturated elements in buf.
	int n_elems = (bitpos + n_bits) / word_size - curr_byte_pos;
	bitpos = (bitpos + n_bits) % max_bitpos;

	// write fibo-encoded values to output buffer for concentrated elements in buf.
	for (int p = 0; p < n_elems; ++p)
	{
		auto &elem = buf[(curr_byte_pos + p) % bufsize];
		obuf.write((char *)&elem, sizeof(elem));
		elem = 0;
	}
}

/**
 * @brief pack elem into buf starting at bitpos, using exactly b_bits bits.
 * @pre num is representable using `b_bits` bits.
 * @pre buf has at least one element
 * @pre buf has at least two elements if bitpos + b_bits > 64.
 *
 * @param b_bits The number of bits to represent elem.
 * @param elem The number to pack, must be representable with at most `b_bits` bits.
 * @param buf The int64_t array where elem should be packed into.
 * @param bitpos The starting point in bits of where to back elem.
 * @return bool: true iff packing of elem saturates buf[0].
 */
template <typename SRC_T, typename DEST_T>
bool pack_int(
	const uint32_t b_bits,
	SRC_T elem,
	DEST_T *buf,
	uint32_t &bitpos)
{
	constexpr int32_t dest_wordsize = sizeof(DEST_T) * 8;
	int32_t shift = dest_wordsize - bitpos - b_bits;
	DEST_T carryover = 0;
	constexpr DEST_T ONE{1};

	if (shift < 0)
	{
		uint32_t r_shift = (b_bits + shift);
		carryover = elem & ((ONE << -shift) - 1);
		*(buf + 1) = carryover << (dest_wordsize + shift);
		elem >>= (-shift);
	}

	*buf |= (((DEST_T)elem) << std::max(0, shift));

	bitpos = (dest_wordsize - shift) % dest_wordsize;
	return (shift <= 0);
}

/**
 * @brief Encode a block of integers using NewPFD.
 * @pre BUF has a size of `pfd_block`.size() / 64 * `b_bits` elements.
 *
 * @param pfd_block The numbers to encode.
 * @param index_gaps Output: delta encoded indices of exceptions in `pfd_block`.
 * @param exceptions Output: The most significant bits of each exception.
 * @param b_bits The number of bits to use per int for the block.
 * @param min_element The smallest element in `pfd_block`.
 * @param BUF The buffer to pack the elements of pfd_block into.
 */
void encode_pfd_block(
	std::vector<int32_t> &pfd_block,
	std::vector<int32_t> &index_gaps,
	std::vector<int32_t> &exceptions,
	const uint32_t b_bits,
	const int32_t min_element,
	PFD_t *BUF)
{
	index_gaps.clear();
	exceptions.clear();

	uint32_t bitpos = 0;
	uint32_t max_elem_bit_mask = (1 << b_bits) - 1;
	uint32_t idx = 0,
			 last_ex_idx = 0;

	bool do_increment_pointer = 0;

	// Store the elements in the primary block in pfd_buf using `b_bits` bits each.
	for (auto &elem : pfd_block)
	{
		elem -= min_element;
		if (elem > max_elem_bit_mask)
		{
			// store the overflowing, most significant bits as exceptions.
			PFD_t exception_bits = (PFD_t)elem >> b_bits;

			exceptions.push_back(elem >> b_bits);
			index_gaps.push_back(idx - last_ex_idx);
			last_ex_idx = idx;

			// store the least b significant bits in the frame.
			elem &= max_elem_bit_mask;
		}

		// pack_int may need more iterations if we do not have PFD_t as 64 bits.
		do_increment_pointer = pack_int<uint32_t>(b_bits, elem, BUF, bitpos);
		BUF += do_increment_pointer;
		++idx;
	}
}

/**
 * @brief Encode a block of `block_size` elements in `pfd_block` and write to of.
 *
 * @param block_size The number of elements in the block.
 * @param pfd_block The elements to encode.
 * @param index_gaps Vector to store delta-encoded indices of exceptions.
 * @param exceptions Vector to store the most significant bits of exceptions.
 * @param fibonacci_buf array used for storing fibonacci encoding.
 * @param b_bits The number of bits each num in `pfd_block` is represented with.
 * @param min_element The smallest element in `pfd_block`.
 * @param of The ostream to write out the encoding.
 */
void new_pfd(
	const size_t block_size,
	std::vector<int32_t> &pfd_block,
	std::vector<int32_t> &index_gaps,
	std::vector<int32_t> &exceptions,
	PFD_t *fibonacci_buf,
	const size_t b_bits,
	const int32_t min_element,
	std::ostream &of)
{
	fibonacci_buf[0] = 0;
	fibonacci_buf[1] = 0;
	fibonacci_buf[2] = 0;

	constexpr size_t wordsize = sizeof(PFD_t) * 8;

	size_t buf_size = block_size / wordsize * b_bits;
	PFD_t *PFD_buf = new PFD_t[buf_size];
	std::fill(PFD_buf, PFD_buf + buf_size, 0ULL);

	PFD_t *const PFD_buf_begin = PFD_buf;
	uint32_t bitpos{0};

	encode_pfd_block(pfd_block, index_gaps, exceptions, b_bits, min_element, PFD_buf);

	size_t n_exceptions = index_gaps.size();

	// For more compact fibonacci encoding, we pack the fibo encoded values together
	fiboEncode<3>(b_bits + 1, fibonacci_buf, bitpos, of);
	fiboEncode<3>(min_element + 1, fibonacci_buf, bitpos, of);
	fiboEncode<3>(n_exceptions + 1, fibonacci_buf, bitpos, of);

	for (const auto &el : index_gaps)
	{
		// we must increment by one, since the first index gap can be zero
		fiboEncode<3>(el + 1, fibonacci_buf, bitpos, of);
	}
	for (const auto &el : exceptions)
	{
		// These are always > 0 since they contain the most significant bits of exception
		fiboEncode<3>(el, fibonacci_buf, bitpos, of);
	}

	flush_fibonacci(fibonacci_buf, bitpos, of);

	of.write((char *)PFD_buf, buf_size * sizeof(PFD_t));
	delete[] PFD_buf;
}

/**
 * @brief Compute the smallest element, and the number of bits required to encode at least 90% of
 * the numbers in pfd_scratch.
 * @pre pfd_scratch contains the numbers to encode in a block
 * @post pfd_scratch may have been reordered.
 * @param block_size The number of elements in pfd_scratch, i.e. the block size.
 * @param pfd_scratch The numbers to encode in a single block using NewPFD.
 * @param min_element Output: The smallest element in `pfd_scratch`.
 * @param b_bits Output: The minimum number of bits that are enough to encode at least ~90% of the elements in `pfd_scratch`.
 */
void compute_pfd_params(
	const size_t block_size,
	std::vector<int32_t> &pfd_scratch,
	int32_t &min_element,
	uint32_t &b_bits)
{
	const size_t nth_elem_idx = (size_t)((double)block_size * 0.9f);
	std::nth_element(pfd_scratch.begin(), pfd_scratch.begin(), pfd_scratch.end());
	min_element = pfd_scratch[0];

	std::nth_element(
		pfd_scratch.begin(),
		pfd_scratch.begin() + nth_elem_idx,
		pfd_scratch.end());

	int32_t nth_max_element = pfd_scratch[nth_elem_idx];
	b_bits = 31 - __builtin_clrsb(nth_max_element - min_element);
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

	constexpr size_t fibonacci_bufsize{cache_size / sizeof(FIBO_t)};
	FIBO_t fibonacci_buf[fibonacci_bufsize];

	std::fill(fibonacci_buf, fibonacci_buf + fibonacci_bufsize, 0);
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
				fiboEncode<fibonacci_bufsize>(1ULL, fibonacci_buf, bit_pos, of);
				fiboEncode<fibonacci_bufsize>(runlen, fibonacci_buf, bit_pos, of);
				runlen = 0;
			}
			fiboEncode<fibonacci_bufsize>(barcode + 1, fibonacci_buf, bit_pos, of);
		}
		last_bc = rows[i].barcode;
	}

	// Take care of the last run of zeros in the delta-encoded barcodes
	if (runlen) {
		fiboEncode<fibonacci_bufsize>(1ULL, fibonacci_buf, bit_pos, of);
		fiboEncode<fibonacci_bufsize>(runlen, fibonacci_buf, bit_pos, of);
	}

	// Write out the last fibonacci element if applicable.
	flush_fibonacci(fibonacci_buf, bit_pos, of);
}

/**
 * @brief Compress UMIs of rows using periodic_delta-runlen(0)-fibonacci encoding and write to `of`.
 *
 * @param rows BUSData array, contains at least `row_count` elements
 * @param row_count The number of UMIs to compress.
 * @param of The ostream for writing the encoding to.
 */
void lossless_compress_umis(BUSData const *const rows, const int row_count, std::ostream &of)
{
	uint64_t last_bc = rows[0].barcode + 1,
			 last_umi = 0,
			 bc, umi, diff;

	constexpr size_t fibonacci_bufsize{cache_size / sizeof(FIBO_t)};
	FIBO_t fibonacci_buf[fibonacci_bufsize];
	std::fill(fibonacci_buf, fibonacci_buf + fibonacci_bufsize, 0U);

	uint32_t bitpos{0};
	uint64_t runlen{0};
	const uint32_t RLE_val{0ULL};

	for (int i = 0; i < row_count; ++i) {
		bc = rows[i].barcode;

		// We must increment umi, since a UMI==0 will confuse the runlength decoder.
		umi = rows[i].UMI + 1;

		if (last_bc != bc) {
			last_umi = 0;
		}

		diff = umi - last_umi;

		if (diff == RLE_val) {
			++runlen;
		}
		else
		{
			// Increment values for fibonacci encoding.
			if (runlen) {
				fiboEncode<fibonacci_bufsize>(RLE_val + 1, fibonacci_buf, bitpos, of);
				fiboEncode<fibonacci_bufsize>(runlen, fibonacci_buf, bitpos, of);
				runlen = 0;
			}
			fiboEncode<fibonacci_bufsize>(diff + 1, fibonacci_buf, bitpos, of);
		}

		last_umi = umi;
		last_bc = bc;
	}

	// Take care of the last run of zeros.
	if (runlen) {
		fiboEncode<fibonacci_bufsize>(RLE_val + 1, fibonacci_buf, bitpos, of);
		fiboEncode<fibonacci_bufsize>(runlen, fibonacci_buf, bitpos, of);
	}

	// Write last bytes if the fibonacci_buf has not been saturated.
	flush_fibonacci(fibonacci_buf, bitpos, of);
}

void lossy_compress_umis(BUSData const *const rows, const int row_count, std::ostream &of)
{
	std::cerr << "BEWARE: Lossy compression\n";

}

/**
 * @brief Compress ECs of rows using NewPFD-fibonacci encoding and write to `of`.
 *
 * @param rows BUSData array, contains at least `row_count` elements
 * @param row_count The number of ECs to compress.
 * @param of The ostream for writing the encoding to.
 */
void compress_ecs(BUSData const *const rows, const int row_count, std::ostream &of)
{
	size_t BLOCK_SIZE{512};
	std::vector<int32_t> index_gaps,
		pfd_scratch,
		pfd_block,
		exceptions;

	index_gaps.reserve(BLOCK_SIZE);
	pfd_scratch.reserve(BLOCK_SIZE);
	pfd_block.reserve(BLOCK_SIZE);
	const size_t fibonacci_bufsize{3};
	PFD_t fibonacci_buf[fibonacci_bufsize];

	int row_index{0};
	int pfd_row_index{0};

	while (row_index < row_count)
	{
		std::fill(fibonacci_buf, fibonacci_buf + fibonacci_bufsize, 0);

		pfd_row_index = 0;
		pfd_scratch.clear();
		pfd_block.clear();

		while (pfd_row_index < BLOCK_SIZE && row_index < row_count)
		{
			pfd_scratch.push_back(rows[row_index].ec);
			pfd_block.push_back(rows[row_index].ec);

			++pfd_row_index;
			++row_index;
		}

		uint32_t b_bits = 0;
		int32_t min_element = 0;
		compute_pfd_params(BLOCK_SIZE, pfd_scratch, min_element, b_bits);
		new_pfd(BLOCK_SIZE, pfd_block, index_gaps, exceptions, fibonacci_buf, b_bits, min_element, of);
	}

	pfd_block.clear();

	// pfd_row_index is then the number of ints in the last block.
	// We signal the end of chunk by encoding b_bits = 0, with the number of elements in the last block as min_element.
	// Since pfd_block is cleared, no additional values will be written (except n_exceptions which adds 2 bits).
	new_pfd(BLOCK_SIZE, pfd_block, index_gaps, exceptions, fibonacci_buf, 0, pfd_row_index, of);
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

	constexpr size_t fibonacci_bufsize{cache_size / sizeof(FIBO_t)};
	FIBO_t fibonacci_buf[fibonacci_bufsize];
	std::fill(fibonacci_buf, fibonacci_buf+fibonacci_bufsize, 0ULL);

	for (int i = 0; i < row_count; ++i)
	{
		count = rows[i].count;

		if (count == RLE_val)
		{
			++runlen;
		} else {
			if (runlen) {
				// Runlength-encode 1s.
				fiboEncode<fibonacci_bufsize>(RLE_val, fibonacci_buf, bitpos, of);
				fiboEncode<fibonacci_bufsize>(runlen, fibonacci_buf, bitpos, of);
				runlen = 0;
			}

			fiboEncode<fibonacci_bufsize>(count, fibonacci_buf, bitpos, of);
		}
	}
	if (runlen)
	{
		// Runlength-encode last run of 1s.
		fiboEncode<fibonacci_bufsize>(RLE_val, fibonacci_buf, bitpos, of);
		fiboEncode<fibonacci_bufsize>(runlen, fibonacci_buf, bitpos, of);
	}

	// Write last bytes when if the fibonacci_buf has not been saturated.
	flush_fibonacci(fibonacci_buf, bitpos, of);
}

/**
 * @brief Compress flags of rows using runlength(0)-fibonacci encoding and write to `of`.
 *
 * @param rows BUSData array, contains at least `row_count` elements
 * @param row_count The number of counts to compress.
 * @param of The ostream for writing the encoding to.
 */
void compress_flags(BUSData const *const rows, const int row_count, std::ostream &of)
{
	const uint32_t RLE_val{0UL};
	uint32_t flag,
		runlen{0},
		bit_pos{0};
	constexpr size_t bufsize{cache_size / sizeof(FIBO_t)};
	FIBO_t buf[bufsize];
	std::fill(buf, buf + bufsize, 0);

	for (int i = 0; i < row_count; ++i)
	{
		flag = rows[i].flags;

		if (flag == RLE_val) {
			++runlen;
		} else {
			// Increment values as fibo cannot encode 0
			if (runlen) {
				// Runlength-encode 0s (incremented).
				fiboEncode<bufsize>(RLE_val + 1, buf, bit_pos, of);
				fiboEncode<bufsize>(runlen, buf, bit_pos, of);
				runlen = 0;
			}
			fiboEncode<bufsize>(flag + 1, buf, bit_pos, of);
		}
	}

	if (runlen) {
		// Runlength-encode last run of 0s (incremented).
		fiboEncode<bufsize>(RLE_val + 1, buf, bit_pos, of);
		fiboEncode<bufsize>(runlen, buf, bit_pos, of);
	}

	// Write last bytes when if the fibonacci_buf has not been saturated.
	flush_fibonacci(buf, bit_pos, of);
}

template <int comp_level>
void compress_barcode_zlib(BUSData const *const rows, const int row_count, std::ostream &of)
{
	uint64_t *buf = new uint64_t[row_count];
	for (int i = 0; i < row_count; ++i)
	{
		buf[i] = rows[i].barcode;
	}
	uLongf src_len = row_count * sizeof(uint64_t);
	uLongf dest_len = row_count * sizeof(uint64_t);
	Bytef *dest = new Bytef[dest_len];
	int status = compress2(dest, &dest_len, (Bytef *)buf, src_len, comp_level);
	if (status == Z_OK)
		of.write((char *)dest, dest_len * sizeof(Bytef));
	else
	{
		std::cerr << "zlib err: " << status << '\n';
	}
}

template <int comp_level>
void compress_UMI_zlib(BUSData const *const rows, const int row_count, std::ostream &of)
{
	uint64_t *buf = new uint64_t[row_count];
	for (int i = 0; i < row_count; ++i)
	{
		buf[i] = rows[i].UMI;
	}
	uLongf src_len = row_count * sizeof(uint64_t);
	uLongf dest_len = row_count * sizeof(uint64_t);
	Bytef *dest = new Bytef[dest_len];
	int status = compress2(dest, &dest_len, (Bytef *)buf, src_len, comp_level);
	if (status == Z_OK)
		of.write((char *)dest, dest_len * sizeof(Bytef));
	else
	{
		std::cerr << "zlib err: " << status << '\n';
	}
}

template <int comp_level>
void compress_EC_zlib(BUSData const * const rows, const int row_count, std::ostream &of)
{
	int32_t *buf = new int32_t[row_count];
	for (int i = 0; i < row_count; ++i){
		buf[i] = rows[i].ec;
	}
	uLongf src_len = row_count * sizeof(int32_t);
	uLongf dest_len = row_count * sizeof(int32_t);
	Bytef *dest = new Bytef[dest_len];
	int status = compress2(dest, &dest_len, (Bytef *)buf, src_len, comp_level);
	if (status == Z_OK)
		of.write((char *)dest, dest_len * sizeof(Bytef));
	else
	{
		std::cerr << "zlib err: " << status << '\n';
	}
}

template<int comp_level>
void compress_count_zlib(BUSData const * const rows, const int row_count, std::ostream &of)
{
	uint32_t *buf = new uint32_t[row_count];
	for (int i = 0; i < row_count; ++i){
		buf[i] = rows[i].count;
	}
	uLongf src_len = row_count * sizeof(uint32_t);
	uLongf dest_len = row_count * sizeof(uint32_t);
	Bytef *dest = new Bytef[dest_len];
	int status = compress2(dest, &dest_len, (Bytef *)buf, src_len, comp_level);
	if (status == Z_OK)
		of.write((char *)dest, dest_len * sizeof(Bytef));
	else
	{
		std::cerr << "zlib err: " << status << '\n';
	}
}
template<int comp_level>
void compress_flags_zlib(BUSData const * const rows, const int row_count, std::ostream &of)
{
	uint32_t *buf = new uint32_t[row_count];
	for (int i = 0; i < row_count; ++i){
		buf[i] = rows[i].flags;
	}
	uLongf src_len = row_count * sizeof(uint32_t);
	uLongf dest_len = row_count * sizeof(uint32_t);
	Bytef *dest = new Bytef[dest_len];
	int status = compress2(dest, &dest_len, (Bytef *)buf, src_len, comp_level);
	if (status == Z_OK)
		of.write((char *)dest, dest_len * sizeof(Bytef));
	else
	{
		std::cerr << "zlib err: " << status << '\n';
	}
}

template <typename T>
void compress_barcode_fibo(BUSData const *const rows, const int row_count, std::ostream &of)
{
	constexpr size_t fibonacci_bufsize(cache_size / sizeof(T));
	T fibonacci_buf[fibonacci_bufsize];
	std::fill(fibonacci_buf, fibonacci_buf + fibonacci_bufsize, 0U);
	uint32_t bitpos{0};
	for (int i = 0; i < row_count; ++i)
	{
		fiboEncode<fibonacci_bufsize>(rows[i].barcode + 1, fibonacci_buf, bitpos, of);
	}
	flush_fibonacci(fibonacci_buf, bitpos, of);
}
template <typename T>
void compress_UMI_fibo(BUSData const *const rows, const int row_count, std::ostream &of)
{
	constexpr size_t fibonacci_bufsize(cache_size / sizeof(T));
	T fibonacci_buf[fibonacci_bufsize];
	std::fill(fibonacci_buf, fibonacci_buf + fibonacci_bufsize, 0U);
	uint32_t bitpos{0};
	for (int i = 0; i < row_count; ++i)
	{
		fiboEncode<fibonacci_bufsize>(rows[i].UMI + 1, fibonacci_buf, bitpos, of);
	}
	flush_fibonacci(fibonacci_buf, bitpos, of);
}
template <typename T>
void compress_EC_fibo(BUSData const *const rows, const int row_count, std::ostream &of)
{
	constexpr size_t fibonacci_bufsize(cache_size / sizeof(T));
	T fibonacci_buf[fibonacci_bufsize];
	std::fill(fibonacci_buf, fibonacci_buf + fibonacci_bufsize, 0U);
	uint32_t bitpos{0};
	for (int i = 0; i < row_count; ++i)
	{
		fiboEncode<fibonacci_bufsize>(rows[i].ec + 1, fibonacci_buf, bitpos, of);
	}
	flush_fibonacci(fibonacci_buf, bitpos, of);
}
template <typename T>
void compress_count_fibo(BUSData const *const rows, const int row_count, std::ostream &of)
{
	constexpr size_t fibonacci_bufsize(cache_size / sizeof(T));
	T fibonacci_buf[fibonacci_bufsize];
	std::fill(fibonacci_buf, fibonacci_buf + fibonacci_bufsize, 0U);
	uint32_t bitpos{0};
	for (int i = 0; i < row_count; ++i)
	{
		fiboEncode<fibonacci_bufsize>(rows[i].count, fibonacci_buf, bitpos, of);
	}
	flush_fibonacci(fibonacci_buf, bitpos, of);
}
template <typename T>
void compress_flags_fibo(BUSData const *const rows, const int row_count, std::ostream &of)
{
	constexpr size_t fibonacci_bufsize(cache_size / sizeof(T));
	T fibonacci_buf[fibonacci_bufsize];
	std::fill(fibonacci_buf, fibonacci_buf + fibonacci_bufsize, 0U);
	uint32_t bitpos{0};
	for (int i = 0; i < row_count; ++i)
	{
		fiboEncode<fibonacci_bufsize>(rows[i].flags + 1, fibonacci_buf, bitpos, of);
	}
	flush_fibonacci(fibonacci_buf, bitpos, of);
}

typedef void (*compress_ptr)(BUSData const *, const int, std::ostream &);

compress_ptr select_zlib_compressor(int col, int lvl)
{
	switch (lvl)
	{
	case 1:
	{
		compress_ptr zlibs[5]{
			&compress_barcode_zlib<1>,
			&compress_UMI_zlib<1>,
			&compress_EC_zlib<1>,
			&compress_count_zlib<1>,
			&compress_flags_zlib<1>};
		return zlibs[col];
	}
	case 2:
	{
		compress_ptr zlibs[5]{
			&compress_barcode_zlib<2>,
			&compress_UMI_zlib<2>,
			&compress_EC_zlib<2>,
			&compress_count_zlib<2>,
			&compress_flags_zlib<2>};
		return zlibs[col];
	}
	case 3:
	{
		compress_ptr zlibs[5]{
			&compress_barcode_zlib<3>,
			&compress_UMI_zlib<3>,
			&compress_EC_zlib<3>,
			&compress_count_zlib<3>,
			&compress_flags_zlib<3>};
		return zlibs[col];
	}
	case 4:
	{
		compress_ptr zlibs[5]{
			&compress_barcode_zlib<4>,
			&compress_UMI_zlib<4>,
			&compress_EC_zlib<4>,
			&compress_count_zlib<4>,
			&compress_flags_zlib<4>};
		return zlibs[col];
	}
	case 5:
	{
		compress_ptr zlibs[5]{
			&compress_barcode_zlib<5>,
			&compress_UMI_zlib<5>,
			&compress_EC_zlib<5>,
			&compress_count_zlib<5>,
			&compress_flags_zlib<5>};
		return zlibs[col];
	}
	case 6:
	{
		compress_ptr zlibs[5]{
			&compress_barcode_zlib<6>,
			&compress_UMI_zlib<6>,
			&compress_EC_zlib<6>,
			&compress_count_zlib<6>,
			&compress_flags_zlib<6>};
		return zlibs[col];
	}
	case 7:
	{
		compress_ptr zlibs[5]{
			&compress_barcode_zlib<7>,
			&compress_UMI_zlib<7>,
			&compress_EC_zlib<7>,
			&compress_count_zlib<7>,
			&compress_flags_zlib<7>};
		return zlibs[col];
	}
	case 8:
	{
		compress_ptr zlibs[5]{
			&compress_barcode_zlib<8>,
			&compress_UMI_zlib<8>,
			&compress_EC_zlib<8>,
			&compress_count_zlib<8>,
			&compress_flags_zlib<8>};
		return zlibs[col];
	}
	case 9:
	{
		compress_ptr zlibs[5]{
			&compress_barcode_zlib<9>,
			&compress_UMI_zlib<9>,
			&compress_EC_zlib<9>,
			&compress_count_zlib<9>,
			&compress_flags_zlib<9>};
		return zlibs[col];
	}
	default:
		return nullptr;
	}
}

uint32_t get_n_rows(std::istream &inf){
	auto header_end = inf.tellg();
	inf.seekg(0, std::ios::end);
	auto file_end = inf.tellg();
	inf.seekg(header_end, std::ios::beg);
	return (file_end - header_end) / 32;
}

void bustools_compress(const Bustools_opt &opt)
{
	BUSHeader h;
	const size_t ROW_SIZE = sizeof(BUSData);
	size_t N = opt.max_memory / ROW_SIZE;
	const size_t chunk_size = (N < opt.chunk_size) ? N : opt.chunk_size;

	void (*funcs[])(BUSData const *, const int, std::ostream &) = {&lossless_compress_umis, &lossy_compress_umis};
	if (opt.lossy_umi)
	{
		std::cerr << "Lossy UMI compression has not been implemented, using lossless instead." << std::endl;
	}

	BUSData *p = new BUSData[chunk_size];

	std::ofstream of;
	std::streambuf *buf = nullptr;
	std::streambuf *headerBuf = nullptr;

	if (opt.stream_out)
	{
		buf = std::cout.rdbuf();
		of.open(opt.temp_files);
		headerBuf = of.rdbuf();
	}
	else
	{
		of.open(opt.output);
		buf = of.rdbuf();
		headerBuf = buf;
	}

	std::ostream outf(buf);
	std::ostream outHeader(headerBuf);

	bool fibonacci_select[5] = {
		((opt.fibo_compress >> 4) & 1) > 0,
		((opt.fibo_compress >> 3) & 1) > 0,
		((opt.fibo_compress >> 2) & 1) > 0,
		((opt.fibo_compress >> 1) & 1) > 0,
		((opt.fibo_compress >> 0) & 1) > 0,
	};

	compress_ptr compressors[5]{
		fibonacci_select[0] ? &compress_barcode_fibo<FIBO_t> : &compress_barcodes,
		fibonacci_select[1] ? &compress_UMI_fibo<FIBO_t> : funcs[false && opt.lossy_umi],
		fibonacci_select[2] ? &compress_EC_fibo<FIBO_t> : &compress_ecs,
		fibonacci_select[3] ? &compress_count_fibo<FIBO_t> : &compress_counts,
		fibonacci_select[4] ? &compress_flags_fibo<FIBO_t> : &compress_flags,
	};

	uint32_t zlib_select = 0;
	for (int i = 0; i < 5; ++i)
	{
		compress_ptr compressor = select_zlib_compressor(i, opt.z_levels[i]);
		if (compressor != nullptr)
		{
			compressors[i] = compressor;
			zlib_select |= (1 << (4-i));
		}
	}

	// TODO: should we really allow multiple files?
	for (const auto &infn : opt.files)
	{
		std::streambuf *inbuf;
		std::ifstream inf;
		uint32_t n_rows = 0;

		if (opt.stream_in)
		{
			inbuf = std::cin.rdbuf();
		}
		else
		{
			inf.open(infn.c_str(), std::ios::binary);
			inbuf = inf.rdbuf();
		}
		std::istream in(inbuf);

		parseHeader(in, h);

		compressed_BUSHeader comp_h;
		comp_h.chunk_size = chunk_size;
		comp_h.lossy_umi = false && opt.lossy_umi;
		comp_h.fibo_zlib_compress = (opt.fibo_compress & ((1 << 5) - 1)) << 5;
		comp_h.fibo_zlib_compress |= zlib_select;

		comp_h.extra_header.text = h.text;
		comp_h.extra_header.version = h.version;
		comp_h.extra_header.bclen = h.bclen;
		comp_h.extra_header.umilen = h.umilen;
		writeCompressedHeader(outHeader, comp_h);

		auto headerLen = outHeader.tellp();
		if (!opt.stream_in)
		{
			n_rows = get_n_rows(in);
			size_t n_blocks = (n_rows - 1) / chunk_size + 1;
			size_t last_chunk_size = n_rows % chunk_size ?: chunk_size;
			uint32_t *block_sizes = new uint32_t[n_blocks];
			outHeader.write((char *)block_sizes, sizeof(uint32_t) * n_blocks);
			delete[] block_sizes;
		}

		uint64_t block_counter = 0;
		size_t last_row_count = 0;

		std::vector<uint32_t> chunk_sizes;
		uint32_t chunk_start = outf.tellp(), chunk_end;

		while (in.good())
		{
			in.read((char *)p, chunk_size * ROW_SIZE);

			size_t row_count = in.gcount() / ROW_SIZE;
			last_row_count = row_count;
			++block_counter;

			for (int i_col = 0; i_col < 5; ++i_col)
			{
				compressors[i_col](p, row_count, outf);
			}
			chunk_end = outf.tellp();
			chunk_sizes.push_back(chunk_end - chunk_start);
			chunk_start = chunk_end;
		}

		outHeader.seekp(0, std::ios_base::beg);

		// Update the compressed header post compression.
		comp_h.last_chunk = last_row_count;
		comp_h.n_chunks = block_counter - 1;

		writeCompressedHeader(outHeader, comp_h);
		outHeader.write((char *)&chunk_sizes[0], sizeof(chunk_sizes[0]) * block_counter);
	}
	delete[] p;
}