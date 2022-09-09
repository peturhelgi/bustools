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

size_t pfd_blocksize = 1024;

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
 * @param bufsize the size of the output buffer.
 * @param buf array of `bufsize` uint64_t. num is fibonacci-encoded in buf.
 * @param bitpos the bit position in buf where the next fibo encoding of num starts.
 * @param obuf the ostream buffer to write to.
 * @return bool true iff encoding the current number does not write outside of buf
 */

const auto fibo_begin = fibo64.begin();
const auto fibo_end = fibo64.end();

template <typename BUF_t>
bool fiboEncode(const uint64_t num, const size_t bufsize, BUF_t *buf, size_t &bitpos)
{
	constexpr uint32_t word_size = sizeof(BUF_t) * 8;

	const size_t max_bitpos = bufsize * word_size;
	constexpr BUF_t ONE{1};

	const uint32_t curr_byte_pos = bitpos / word_size;

	uint64_t remainder = num;

	// the ith fibonacci number is the largest fibo not greater than remainder
	auto i = std::upper_bound(fibo_begin, fibo_end, remainder) - 1;

	const uint32_t n_bits = (i - fibo_begin) + 2;

	// Encoding the current number would write out of bounds of buf.
	if (bitpos + n_bits > max_bitpos)
		return false;

	uint32_t next_bit_pos = bitpos + n_bits - 1,
			 bit_offset = next_bit_pos % word_size,
			 buf_offset = (next_bit_pos / word_size) % bufsize;

	// Set the stop bit.
	// buf[buf_offset] |= shifted_64[word_size - 1 - bit_offset];
	buf[buf_offset] |= ONE << (word_size - 1 - bit_offset);

	++i;
	while (remainder > 0)
	{
		i = std::upper_bound(fibo_begin, i+1, remainder) - 1;
		next_bit_pos = bitpos + (i - fibo_begin);
		buf_offset = (next_bit_pos / word_size) % bufsize;
		bit_offset = next_bit_pos % word_size;

		// buf[buf_offset] |= shifted_64[word_size - 1 - bit_offset];
		buf[buf_offset] |= ONE << (word_size - 1 - bit_offset);
		remainder -= *i;
	}
	bitpos += n_bits;
	return true;
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
		// TODO: is it better to use pack_num?
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
size_t new_pfd(
	const size_t block_size,
	std::vector<int32_t> &pfd_block,
	std::vector<int32_t> &index_gaps,
	std::vector<int32_t> &exceptions,
	PFD_t *fibonacci_buf,
	const size_t b_bits,
	const int32_t min_element,
	size_t fibonacci_bufsize,
	PFD_t *PFD_buf
)
{
	bool success = true;
	constexpr size_t wordsize = sizeof(PFD_t) * 8;

	size_t buf_size = block_size * b_bits / wordsize;
	std::fill(PFD_buf, PFD_buf + buf_size, 0ULL);

	size_t bitpos{0};

	encode_pfd_block(pfd_block, index_gaps, exceptions, b_bits, min_element, PFD_buf);

	size_t n_exceptions = index_gaps.size();
	// For more compact fibonacci encoding, we pack the fibo encoded values together
	success &= fiboEncode(b_bits + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
	success &= fiboEncode(min_element + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
	success &= fiboEncode(n_exceptions + 1, fibonacci_bufsize, fibonacci_buf, bitpos);

	for (const auto &el : index_gaps)
	{
		// we must increment by one, since the first index gap can be zero
		success &= fiboEncode(el + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
	}
	for (const auto &el : exceptions)
	{
		// These are always > 0 since they contain the most significant bits of exception
		success &= fiboEncode(el, fibonacci_bufsize, fibonacci_buf, bitpos);
	}

	size_t n_elems = bitpos / wordsize + (bitpos % wordsize > 0);
	success &= (n_elems + buf_size <= fibonacci_bufsize);

	std::memcpy(fibonacci_buf + n_elems, PFD_buf, buf_size * sizeof(PFD_t) * success);
	n_elems += buf_size;

	return n_elems * success;
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
	const size_t nth_elem_idx = (size_t)((double)block_size * 0.9);
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
 * @return bool true iff encoding does not go out of bounds of obuf
 */
bool compress_barcodes(BUSData const *const rows, const int row_count, char *obuf, const size_t &obuf_size, size_t &global_bufpos)
{
	bool success = true;
	uint64_t barcode = 0,
			 last_bc = 0,
			 runlen = 0;
	size_t wordsize = sizeof(FIBO_t) * 8;

	const size_t fibonacci_bufsize = (obuf_size - global_bufpos) / sizeof(FIBO_t);
	FIBO_t *fibonacci_buf = (FIBO_t *)(obuf + global_bufpos);
	size_t bitpos{0};

	for (int i = 0; i < row_count && success; ++i)
	{
		barcode = rows[i].barcode;

		// delta encoding
		barcode -= last_bc;

		// Runlength encoding of zeros
		if (barcode == 0)
		{
			++runlen;
		}
		else
		{
			// Increment values as fibo cannot encode 0
			if (runlen)
			{
				success &= fiboEncode(1ULL, fibonacci_bufsize, fibonacci_buf, bitpos);
				success &= fiboEncode(runlen, fibonacci_bufsize, fibonacci_buf, bitpos);
				runlen = 0;
			}
			success &= fiboEncode(barcode + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
		}
		last_bc = rows[i].barcode;
	}

	// Take care of the last run of zeros in the delta-encoded barcodes
	if (runlen)
	{
		success &= fiboEncode(1ULL, fibonacci_bufsize, fibonacci_buf, bitpos);
		success &= fiboEncode(runlen, fibonacci_bufsize, fibonacci_buf, bitpos);
	}

	global_bufpos += (bitpos / wordsize + (bitpos % wordsize > 0)) * sizeof(FIBO_t);
	return success;
}

/**
 * @brief Compress UMIs of rows using periodic_delta-runlen(0)-fibonacci encoding and write to `of`.
 *
 * @param rows BUSData array, contains at least `row_count` elements
 * @param row_count The number of UMIs to compress.
 * @param of The ostream for writing the encoding to.
 * @return bool true iff encoding does not go out of bounds of obuf
 */
bool lossless_compress_umis(BUSData const *const rows, const int row_count, char *obuf, const size_t &obuf_size, size_t &global_bufpos)
{
	bool success = true;
	uint64_t last_bc = rows[0].barcode + 1,
			 last_umi = 0,
			 bc, umi, diff;

	size_t wordsize = sizeof(FIBO_t) * 8;

	const size_t fibonacci_bufsize = (obuf_size - global_bufpos) / sizeof(FIBO_t);
	FIBO_t *fibonacci_buf = (FIBO_t *)(obuf + global_bufpos);
	size_t bitpos{0};

	uint64_t runlen{0};
	const uint32_t RLE_val{0ULL};

	for (int i = 0; i < row_count && success; ++i)
	{
		bc = rows[i].barcode;

		// We must increment umi, since a UMI==0 will confuse the runlength decoder.
		umi = rows[i].UMI + 1;

		if (last_bc != bc)
		{
			last_umi = 0;
		}

		diff = umi - last_umi;

		if (diff == RLE_val)
		{
			++runlen;
		}
		else
		{
			// Increment values for fibonacci encoding.
			if (runlen)
			{
				success &= fiboEncode(RLE_val + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
				success &= fiboEncode(runlen, fibonacci_bufsize, fibonacci_buf, bitpos);

				runlen = 0;
			}

			success &= fiboEncode(diff + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
		}

		last_umi = umi;
		last_bc = bc;
	}

	// Take care of the last run of zeros.
	if (runlen)
	{
		success &= fiboEncode(RLE_val + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
		success &= fiboEncode(runlen, fibonacci_bufsize, fibonacci_buf, bitpos);
	}

	global_bufpos += (bitpos / wordsize + (bitpos % wordsize > 0)) * sizeof(FIBO_t);

	return success;
}

bool lossy_compress_umis(BUSData const *const rows, const int row_count, char *obuf, const size_t &obuf_size, size_t &global_bufpos)
{
	bool success = true;
	std::cerr << "BEWARE: Lossy compression\n";
	return success;
}

/**
 * @brief Compress ECs of rows using NewPFD-fibonacci encoding and write to `of`.
 *
 * @param rows BUSData array, contains at least `row_count` elements
 * @param row_count The number of ECs to compress.
 * @param of The ostream for writing the encoding to.
 * @return bool true iff encoding does not go out of bounds of obuf
 */
bool compress_ecs(BUSData const *const rows, const int row_count, char *obuf, const size_t &obuf_size, size_t &global_bufpos)
{
	bool success = true;
	size_t BLOCK_SIZE{pfd_blocksize};
	size_t wordsize = sizeof(PFD_t) * 8;
	size_t buf_offset = 0;

	std::vector<int32_t> index_gaps,
		pfd_scratch,
		pfd_block,
		exceptions;

	// todo: We might be able to speed up by creating the primary array here.
	size_t max_size_block = BLOCK_SIZE * sizeof(int32_t) / sizeof(PFD_t);
	PFD_t *primary_block = new PFD_t[max_size_block];

	exceptions.reserve(BLOCK_SIZE);
	index_gaps.reserve(BLOCK_SIZE);
	pfd_scratch.reserve(BLOCK_SIZE);
	pfd_block.reserve(BLOCK_SIZE);

	const size_t fibonacci_bufsize = (obuf_size - global_bufpos) / sizeof(PFD_t);
	PFD_t *fibonacci_buf = (PFD_t *)(obuf + global_bufpos);

	int row_index{0};
	int pfd_row_index{0};
	size_t elems_written = 0;
	size_t byte_count = 0;

	while (row_index < row_count && success)
	{
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

		// we don't want to reset the fibonacci bytes, as this is or primary buffer.
		elems_written = new_pfd(BLOCK_SIZE, pfd_block, index_gaps, exceptions, fibonacci_buf + buf_offset, b_bits, min_element, fibonacci_bufsize - buf_offset, primary_block);

		success &= (elems_written > 0);
		buf_offset += elems_written;
		byte_count += elems_written * sizeof(PFD_t);
	}

	pfd_block.clear();

	// pfd_row_index is then the number of ints in the last block.
	// We signal the end of chunk by encoding b_bits = 0, with the number of elements in the last block as min_element.
	// Since pfd_block is cleared, no additional values will be written (except n_exceptions which adds 2 bits).

	elems_written = new_pfd(BLOCK_SIZE, pfd_block, index_gaps, exceptions, fibonacci_buf + buf_offset, 0, pfd_row_index, fibonacci_bufsize - buf_offset, primary_block);

	success &= (elems_written > 0);
	byte_count += elems_written * sizeof(PFD_t);
	global_bufpos += byte_count;

	delete[] primary_block;
	return success;
}

/**
 * @brief Compress counts of rows using runlength(1)-fibonacci encoding and write to `of`.
 *
 * @param rows BUSData array, contains at least `row_count` elements
 * @param row_count The number of counts to compress.
 * @param of The ostream for writing the encoding to.
 * @return bool true iff encoding does not go out of bounds of obuf
 */
bool compress_counts(BUSData const *const rows, const int row_count, char *obuf, const size_t &obuf_size, size_t &global_bufpos)
{
	bool success = true;
	const uint32_t RLE_val{1UL};
	uint32_t count,
		runlen{0};

	size_t wordsize = sizeof(FIBO_t) * 8;

	const size_t fibonacci_bufsize = (obuf_size - global_bufpos) / sizeof(FIBO_t);
	FIBO_t *fibonacci_buf = (FIBO_t *)(obuf + global_bufpos);
	size_t bitpos{0};

	for (int i = 0; i < row_count && success; ++i)
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
				success &= fiboEncode(RLE_val, fibonacci_bufsize, fibonacci_buf, bitpos);
				success &= fiboEncode(runlen, fibonacci_bufsize, fibonacci_buf, bitpos);
				runlen = 0;
			}
			success &= fiboEncode(count, fibonacci_bufsize, fibonacci_buf, bitpos);
		}
	}
	if (runlen)
	{
		// Runlength-encode last run of 1s.
		success &= fiboEncode(RLE_val, fibonacci_bufsize, fibonacci_buf, bitpos);
		success &= fiboEncode(runlen, fibonacci_bufsize, fibonacci_buf, bitpos);
	}

	global_bufpos += (bitpos / wordsize + (bitpos % wordsize > 0)) * sizeof(FIBO_t);
	return success;
}

/**
 * @brief Compress flags of rows using runlength(0)-fibonacci encoding and write to `of`.
 *
 * @param rows BUSData array, contains at least `row_count` elements
 * @param row_count The number of counts to compress.
 * @param of The ostream for writing the encoding to.
 * @return bool true iff encoding does not go out of bounds of obuf
 */
bool compress_flags(BUSData const *const rows, const int row_count, char *obuf, const size_t &obuf_size, size_t &global_bufpos)
{
	bool success = true;
	const uint32_t RLE_val{0UL};
	uint32_t flag,
		runlen{0};

	size_t wordsize = sizeof(FIBO_t) * 8;

	const size_t fibonacci_bufsize = (obuf_size - global_bufpos) / sizeof(FIBO_t);
	FIBO_t *fibonacci_buf = (FIBO_t *)(obuf + global_bufpos);
	size_t bitpos{0};
	// don't need to fill with zeros as that is done prior.

	for (int i = 0; i < row_count && success; ++i)
	{
		flag = rows[i].flags;

		if (flag == RLE_val)
		{
			++runlen;
		}
		else
		{
			// Increment values as fibo cannot encode 0
			if (runlen)
			{
				// Runlength-encode 0s (incremented).
				success &= fiboEncode(RLE_val + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
				success &= fiboEncode(runlen, fibonacci_bufsize, fibonacci_buf, bitpos);
				runlen = 0;
			}
			success &= fiboEncode(flag + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
		}
	}

	if (runlen)
	{
		// Runlength-encode last run of 0s (incremented).
		success &= fiboEncode(RLE_val + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
		success &= fiboEncode(runlen, fibonacci_bufsize, fibonacci_buf, bitpos);
	}

	global_bufpos += (bitpos / wordsize + (bitpos % wordsize > 0)) * sizeof(FIBO_t);
	return success;
}

template <int comp_level>
bool compress_barcode_zlib(BUSData const *const rows, const int row_count, char *obuf, const size_t &obuf_size, size_t &global_bufpos)
{
	bool success = true;
	typedef uint64_t T;

	T *src_buf = new T[row_count];
	for (int i = 0; i < row_count; ++i)
	{
		src_buf[i] = rows[i].barcode;
	}

	uLongf src_len = row_count * sizeof(T);
	uLongf dest_len = (obuf_size - global_bufpos) / sizeof(Bytef);

	Bytef *dest = (Bytef *)(obuf + global_bufpos);
	int status = compress2(dest, &dest_len, (Bytef *)src_buf, src_len, comp_level);
	if (status != Z_OK)
	{
		success = false;
		std::cerr << "zlib err: " << status << '\n';
	}

	global_bufpos += dest_len;
	delete[] src_buf;
	return success;
}

template <int comp_level>
bool compress_UMI_zlib(BUSData const *const rows, const int row_count, char *obuf, const size_t &obuf_size, size_t &global_bufpos)
{
	bool success = true;
	typedef uint64_t T;

	T *src_buf = new T[row_count];
	for (int i = 0; i < row_count; ++i)
	{
		src_buf[i] = rows[i].UMI;
	}
	uLongf src_len = row_count * sizeof(T);
	uLongf dest_len = (obuf_size - global_bufpos) / sizeof(Bytef);

	Bytef *dest = (Bytef *)(obuf + global_bufpos);
	int status = compress2(dest, &dest_len, (Bytef *)src_buf, src_len, comp_level);
	if (status != Z_OK)
	{
		success = false;
		std::cerr << "zlib err: " << status << '\n';
	}

	global_bufpos += dest_len;
	delete[] src_buf;
	return success;
}

template <int comp_level>
bool compress_EC_zlib(BUSData const *const rows, const int row_count, char *obuf, const size_t &obuf_size, size_t &global_bufpos)
{
	bool success = true;
	typedef int32_t T;

	T *src_buf = new T[row_count];
	for (int i = 0; i < row_count; ++i)
	{
		src_buf[i] = rows[i].ec;
	}
	uLongf src_len = row_count * sizeof(T);
	uLongf dest_len = (obuf_size - global_bufpos) / sizeof(Bytef);

	Bytef *dest = (Bytef *)(obuf + global_bufpos);
	int status = compress2(dest, &dest_len, (Bytef *)src_buf, src_len, comp_level);
	if (status != Z_OK)
	{
		success = false;
		std::cerr << "zlib err: " << status << '\n';
	}

	global_bufpos += dest_len;
	delete[] src_buf;
	return success;
}

template <int comp_level>
bool compress_count_zlib(BUSData const *const rows, const int row_count, char *obuf, const size_t &obuf_size, size_t &global_bufpos)
{
	bool success = true;
	typedef uint32_t T;

	T *src_buf = new T[row_count];
	for (int i = 0; i < row_count; ++i)
	{
		src_buf[i] = rows[i].count;
	}
	uLongf src_len = row_count * sizeof(T);
	uLongf dest_len = (obuf_size - global_bufpos) / sizeof(Bytef);

	Bytef *dest = (Bytef *)(obuf + global_bufpos);
	int status = compress2(dest, &dest_len, (Bytef *)src_buf, src_len, comp_level);
	if (status != Z_OK)
	{
		success = false;
		std::cerr << "zlib err: " << status << '\n';
	}

	global_bufpos += dest_len;
	delete[] src_buf;
	return success;
}
template <int comp_level>
bool compress_flags_zlib(BUSData const *const rows, const int row_count, char *obuf, const size_t &obuf_size, size_t &global_bufpos)
{
	bool success = true;

	typedef uint32_t T;
	T *src_buf = new T[row_count];
	for (int i = 0; i < row_count; ++i)
	{
		src_buf[i] = rows[i].flags;
	}
	uLongf src_len = row_count * sizeof(T);
	uLongf dest_len = (obuf_size - global_bufpos) / sizeof(Bytef);

	Bytef *dest = (Bytef *)(obuf + global_bufpos);
	int status = compress2(dest, &dest_len, (Bytef *)src_buf, src_len, comp_level);
	if (status != Z_OK)
	{
		std::cerr << "zlib err: " << status << '\n';
		success = false;
	}

	global_bufpos += dest_len;
	delete[] src_buf;
	return success;
}

template <typename T>
bool compress_barcode_fibo(BUSData const *const rows, const int row_count, char *obuf, const size_t &obuf_size, size_t &global_bufpos)
{
	bool success = true;
	size_t bitpos{0};
	size_t wordsize = sizeof(T) * 8;

	const size_t fibonacci_bufsize = (obuf_size - global_bufpos) / sizeof(T);
	T *fibonacci_buf = (T *)(obuf + global_bufpos);

	for (int i = 0; i < row_count; ++i)
	{
		success &= fiboEncode(rows[i].barcode + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
	}

	global_bufpos += (bitpos / wordsize + (bitpos % wordsize > 0)) * sizeof(T);
	return success;
}
template <typename T>
bool compress_UMI_fibo(BUSData const *const rows, const int row_count, char *obuf, const size_t &obuf_size, size_t &global_bufpos)
{
	bool success = true;
	size_t bitpos{0};
	size_t wordsize = sizeof(T) * 8;

	const size_t fibonacci_bufsize = (obuf_size - global_bufpos) / sizeof(T);
	T *fibonacci_buf = (T *)(obuf + global_bufpos);

	for (int i = 0; i < row_count; ++i)
	{
		success &= fiboEncode(rows[i].UMI + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
	}

	global_bufpos += (bitpos / wordsize + (bitpos % wordsize > 0)) * sizeof(T);
	return success;
}
template <typename T>
bool compress_EC_fibo(BUSData const *const rows, const int row_count, char *obuf, const size_t &obuf_size, size_t &global_bufpos)
{
	bool success = true;
	size_t bitpos{0};
	size_t wordsize = sizeof(T) * 8;

	const size_t fibonacci_bufsize = (obuf_size - global_bufpos) / sizeof(T);
	T *fibonacci_buf = (T *)(obuf + global_bufpos);

	for (int i = 0; i < row_count; ++i)
	{
		success &= fiboEncode(rows[i].ec + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
	}

	global_bufpos += (bitpos / wordsize + (bitpos % wordsize > 0)) * sizeof(T);
	return success;
}
template <typename T>
bool compress_count_fibo(BUSData const *const rows, const int row_count, char *obuf, const size_t &obuf_size, size_t &global_bufpos)
{
	bool success = true;
	size_t bitpos{0};
	size_t wordsize = sizeof(T) * 8;

	const size_t fibonacci_bufsize = (obuf_size - global_bufpos) / sizeof(T);
	T *fibonacci_buf = (T *)(obuf + global_bufpos);

	for (int i = 0; i < row_count; ++i)
	{
		success &= fiboEncode(rows[i].count, fibonacci_bufsize, fibonacci_buf, bitpos);
	}

	global_bufpos += (bitpos / wordsize + (bitpos % wordsize > 0)) * sizeof(T);
	return success;
}
template <typename T>
bool compress_flags_fibo(BUSData const *const rows, const int row_count, char *obuf, const size_t &obuf_size, size_t &global_bufpos)
{
	bool success = true;
	size_t bitpos{0};
	size_t wordsize = sizeof(T) * 8;

	const size_t fibonacci_bufsize = (obuf_size - global_bufpos) / sizeof(T);
	T *fibonacci_buf = (T *)(obuf + global_bufpos);

	for (int i = 0; i < row_count; ++i)
	{
		success &= fiboEncode(rows[i].flags + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
	}

	global_bufpos += (bitpos / wordsize + (bitpos % wordsize > 0)) * sizeof(T);
	return success;
}

typedef bool (*compress_ptr)(BUSData const *, const int, char *, const size_t &, size_t &);

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

int get_target_file_type(const std::string &filename)
{
	std::ifstream inf(filename);
	if (!inf.is_open())
	{
		return -2;
	}
	char magic[5];
	int ret = -1;
	inf.read((char *)&magic[0], 4);

	magic[4] = '\0';
	if (std::strcmp(&magic[0], "BUS\0") == 0)
	{
		return 0;
	}
	if (std::strcmp(&magic[0], "BUS\1\0") == 0)
	{
		return 1;
	}

	if (std::strcmp(&magic[0], "0\t0\n") == 0)
	{
		return 2;
	}

	if (std::strcmp(&magic[0], "BEC\0") == 0)
	{
		return 3;
	}
	return -1;
}

uint32_t select_compressors(const Bustools_opt &opt, compress_ptr compressors[5]){
	
	compress_ptr fibonacci_compressors[5]{
		&compress_barcode_fibo<FIBO_t>,
		&compress_UMI_fibo<FIBO_t>,
		&compress_EC_fibo<FIBO_t>,
		&compress_count_fibo<FIBO_t>,
		&compress_flags_fibo<FIBO_t>,
	};
	compress_ptr bustools_compressors[5]{
		&compress_barcodes,
		(opt.lossy_umi ? &lossy_compress_umis : &lossless_compress_umis),
		&compress_ecs,
		&compress_counts,
		&compress_flags,
	};

	uint32_t zlib_select = 0;
	for (int i = 0; i < 5; ++i){
		compressors[i] = bustools_compressors[i];
		if((opt.fibo_compress >> (4-i) & 1) > 0)
		{
			compressors[i] = fibonacci_compressors[i];
		}
		else{
			compress_ptr comp = select_zlib_compressor(i, opt.z_levels[i]);
			if(comp != nullptr){
				compressors[i] = comp;
				zlib_select |= (1 << (4 - i));
			}
		}
	}

	uint32_t fibo_select = opt.fibo_compress & ((1 << 5) - 1);
	return (fibo_select << 5) | zlib_select;
}

void compress_busfile(const Bustools_opt &opt, std::ostream &outf, std::ostream &outHeaderf, std::streambuf *inbuf, BUSHeader &h)
{
	constexpr size_t ROW_SIZE = sizeof(BUSData);

	size_t N = opt.max_memory / ROW_SIZE;
	const size_t chunk_size = (N < opt.chunk_size) ? N : opt.chunk_size;

	compress_ptr compressors[5];
	uint32_t compressor_selection = select_compressors(opt, compressors);
	std::istream in(inbuf);
	
	bool parsed = parseHeader(in, h);
	
	compressed_BUSHeader comp_h;
	comp_h.chunk_size = chunk_size;
	comp_h.lossy_umi = opt.lossy_umi;
	comp_h.fibo_zlib_compress = compressor_selection;

	comp_h.extra_header.text = h.text;
	comp_h.extra_header.version = h.version;
	comp_h.extra_header.bclen = h.bclen;
	comp_h.extra_header.umilen = h.umilen;

	pfd_blocksize = opt.pfd_blocksize;
	comp_h.pfd_blocksize = pfd_blocksize;

	writeCompressedHeader(outHeaderf, comp_h);

	std::vector<uint32_t> block_sizes;

	// 6 * chunk_sizes is usually good enough, but we make it a multiple of 8;
	size_t bufsize = (6 * chunk_size / 8) * 8;
	size_t bufpos = 0;
	size_t buf_checkpoint = 0;
	size_t row_count = 0;

	uint64_t block_counter = 0;
	uint64_t block_header = 0;

	try
	{
		BUSData *busdata = new BUSData[chunk_size];
		char *buffer = new char[bufsize];
		std::fill(buffer, buffer + bufsize, 0);
		while (in.good())
		{
			in.read((char *)busdata, chunk_size * ROW_SIZE);
			row_count = in.gcount() / ROW_SIZE;

			for (int i_col = 0; i_col < 5; ++i_col)
			{
				bool success = compressors[i_col](busdata, row_count, buffer, bufsize, bufpos);
				if (!success)
				{
					bufsize *= 2;

					char *newbuf = new char[bufsize];
					std::memcpy(newbuf, buffer, buf_checkpoint);
					std::fill(newbuf + buf_checkpoint, newbuf + bufsize, 0);

					delete[] buffer;
					buffer = newbuf;

					bufpos = buf_checkpoint;
					--i_col;
				}
				buf_checkpoint = bufpos;
			}

			block_header = bufpos << 30;
			block_header |= row_count;

			outf.write((char *)&block_header, sizeof(block_header));
			outf.write(buffer, bufpos);

			std::fill(buffer, buffer + bufpos, 0);

			block_sizes.push_back(bufpos);
			bufpos = 0;
			++block_counter;
		}
		block_header = 0;
		outf.write((char *)&block_header, sizeof(block_header));

		// todo: move these to index file or back of file
		// todo: put first barcodes in index file as well

		outHeaderf.seekp(0, std::ios_base::beg);
		comp_h.last_chunk = row_count;
		comp_h.n_chunks = block_counter - 1;

		writeCompressedHeader(outHeaderf, comp_h);

		delete[] busdata;
		delete[] buffer;
	}
	catch (const std::bad_alloc &ex)
	{
		std::cerr << "Unable to allocate buffer\n"
				  << ex.what() << std::endl;
	}
}
template <typename T>
void pack_ec_row_to_file(
	const std::vector<int32_t> &ecs,
	const size_t bufsize,
	T *buf,
	size_t &bitpos,
	std::ostream &of)
{
	bool success = true;
	// Diffs must be incremented since ec == 0 is valid, although rare.
	constexpr int32_t RL_VAL{2};

	size_t n_elems{ecs.size()};

	uint32_t diff,
		last_ec{0},
		runlen{0};

	success &= fiboEncode<T>(n_elems, bufsize, buf, bitpos);

	for (const auto &ec : ecs)
	{
		diff = ec - last_ec + 1;
		if (diff == RL_VAL)
		{
			++runlen;
		}
		else
		{
			if (runlen)
			{
				success &= fiboEncode<T>(RL_VAL, bufsize, buf, bitpos);
				success &= fiboEncode<T>(runlen, bufsize, buf, bitpos);

				runlen = 0;
			}
			success &= fiboEncode<T>(diff, bufsize, buf, bitpos);
		}
		last_ec = ec;
	}
	if (runlen)
	{
		success &= fiboEncode<T>(RL_VAL, bufsize, buf, bitpos);
		success &= fiboEncode<T>(runlen, bufsize, buf, bitpos);
	}
}

template <typename T = uint16_t>
void compress_ec_matrix(const std::string &filename, BUSHeader &h)
{
	// first m rows map to themselves. In the file, we count how many such rows there are.
	// For the rest, we don't need the ids, only the delta+runlen-1+int-coding of the ECs
	parseECs(filename, h);
	std::cerr << "Done parsing ecs" << std::endl;
	auto &ecs = h.ecs;

	uint32_t lo = 0, hi = ecs.size() - 1, mid;

	// Assume all i->i rows are at the beginning of file.
	while (lo < hi)
	{
		// |  ecs[i] = [i]  |   ?   |  ecs[i] != [i]  |
		//  ^              ^       ^                   ^
		//  0              lo     hi                 ecs.size()
		mid = lo + (hi - lo + 1) / 2;
		if (ecs.at(mid)[0] != mid || ecs.at(mid).size() != 1)
		{
			hi = mid - 1;
		}
		else
		{
			lo = mid;
		}
	}

	uint32_t num_identities = lo + 1;
	std::string filename_out = filename + 'z';

	constexpr size_t fibonacci_bufsize{24 / sizeof(T)};

	T fibonacci_buf[fibonacci_bufsize];
	std::fill(fibonacci_buf, fibonacci_buf + fibonacci_bufsize, 0);

	std::ofstream out(filename_out.c_str(), std::ios::binary);
	std::ostream of(out.rdbuf());

	of.write("BEC\0", 4);
	of.write((char *)&num_identities, sizeof(num_identities));
	num_identities = ecs.size() - num_identities;
	of.write((char *)&num_identities, sizeof(num_identities));

	size_t bitpos{0};
	const auto ecs_end = ecs.end();

	for (auto ecs_it = ecs.begin() + lo + 1; ecs_it < ecs_end; ++ecs_it)
	{
		auto &ecs_list = *ecs_it;
		pack_ec_row_to_file(ecs_list, fibonacci_bufsize, fibonacci_buf, bitpos, of);
	}

	flush_fibonacci<T>(fibonacci_buf, bitpos, of);
}

void bustools_compress(const Bustools_opt &opt)
{
	BUSHeader h;

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

	// TODO: should we really allow multiple files?
	for (const auto &infn : opt.files)
	{
		std::streambuf *inbuf;
		std::ifstream inf;

		if (opt.stream_in)
		{
			inbuf = std::cin.rdbuf();
		}
		else
		{
			int target_file_type = get_target_file_type(infn);

			inf.open(infn.c_str(), std::ios::binary);
			inbuf = inf.rdbuf();

			switch (target_file_type)
			{
			// compress busfile
			case 0:
				std::cerr << "Compressing BUS file " << infn << "\n";
				compress_busfile(opt, outf, outHeader, inbuf, h);
				break;
			// decompress busz file
			case 1:
				std::cerr << "Warning: The file " << infn << " is a compressed BUS file. Skipping.\n";
				break;
			case 2:
				std::cerr << "Compressing matrix file " << infn << '\n';
				compress_ec_matrix(infn, h);
				break;
			// read compressed matrix.ecz
			case 3:
				std::cerr << "Warning: The file " << infn << " is a compressed EC matrix file. Skipping.\n";
				break;
			case -2:
				std::cerr << "Error: Unable to open file " << infn << '\n';
				break;
			default:

				std::cerr << "Warning: Unknown file type. Skipping compression on" << infn << "\n";
				break;
			}
		}
	}
}
