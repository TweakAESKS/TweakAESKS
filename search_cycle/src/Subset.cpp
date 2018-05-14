/*
 * Subset.cpp
 *
 * Copyright 2016 Patrick Derbez <patrick.derbez@irisa.fr>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 *
 *
 */


#include <nmmintrin.h>
#include "Subset.hpp"


using namespace std;


Subset & Subset::operator=(unsigned int x)
{
	auto const bound = (x >> 6);
	if (u_size_array <= bound)
	{
		delete[] ptr_array;
		ptr_array = new uint64_t [bound+1];
	}
	u_size_array = bound+1;
	unsigned int i = 0;
	for (; i < bound; ++i) ptr_array[i] = 0;
	ptr_array[i] = uint64_t (1) << (x & 0x3F);
	return *this;
}

Subset & Subset::operator+=(Subset const & s)
{
	if (u_size_array < s.u_size_array)
	{
		uint64_t * __restrict__ ptr = new uint64_t [s.u_size_array];
		unsigned int i = 0;
		for (; i < u_size_array; ++i) ptr[i] = ptr_array[i] | s.ptr_array[i];
		for (; i < s.u_size_array; ++i) ptr[i] = s.ptr_array[i];
		delete[] ptr_array;
		ptr_array = ptr;
		u_size_array = s.u_size_array;
	}
	else
	{
		for (unsigned int i = 0; i < s.u_size_array; ++i) ptr_array[i] |= s.ptr_array[i];
	}
	return *this;
}


void Subset::emplace(unsigned int x)
{
	auto const bound = (x >> 6);
	if (u_size_array <= bound) {
		uint64_t * __restrict__ ptr = new uint64_t [bound+1];
		unsigned int i = 0;
		for (; i < u_size_array; ++i) ptr[i] = ptr_array[i];
		for (; i < bound; ++i) ptr[i] = 0;
		ptr[i] = uint64_t (1) << (x & 0x3F);
		delete[] ptr_array;
		ptr_array = ptr;
		u_size_array = bound + 1;
	}
	else ptr_array[bound] |= uint64_t (1) << (x & 0x3F);
}




Subset & Subset::operator^=(Subset const & s)
{
	if (u_size_array < s.u_size_array)
	{
		uint64_t * __restrict__ ptr = new uint64_t [s.u_size_array];
		unsigned int i = 0;
		for (; i < u_size_array; ++i) ptr[i] = ptr_array[i] ^ s.ptr_array[i];
		for (; i < s.u_size_array; ++i) ptr[i] = s.ptr_array[i];
		delete[] ptr_array;
		ptr_array = ptr;
		u_size_array = s.u_size_array;
	}
	else
	{
		for (unsigned int i = 0; i < s.u_size_array; ++i) ptr_array[i] ^= s.ptr_array[i];
	}
	return *this;
}

Subset & Subset::operator^=(Subset && s)
{
	if (u_size_array < s.u_size_array) swap(*this,s);
	for (unsigned int i = 0; i < s.u_size_array; ++i) ptr_array[i] ^= s.ptr_array[i];
	return *this;
}

unsigned popcount64d(uint64_t x)
{
    unsigned count;
    for (count=0; x; count++) {
        x &= x - 1;
			}
    return count;
}

unsigned int Subset::size() const
{
	unsigned int res = 0;
	//for (unsigned int i = 0; i < u_size_array; ++i) res += _mm_popcnt_u64(ptr_array[i]);
	for (unsigned int i = 0; i < u_size_array; ++i) res += popcount64d(ptr_array[i]);
	return res;
}

bool Subset::empty() const
{
	for (unsigned i = 0; i < u_size_array; ++i) if (ptr_array[i] != 0) return false;
	return true;
}



std::vector<unsigned int> Subset::getElements() const
{
	std::vector<unsigned int> res;
	for (unsigned int i = 0; i < u_size_array; ++i)
	{
		auto x = ptr_array[i];
		unsigned int const n = 64*i;
		unsigned int m = 0;
		while (x != 0)
		{
			while ((x & 1) == 0) {
				x >>= 1;
				++m;
			}
			res.emplace_back(n + m);
			x >>= 1;
			++m;
		}
	}
	return res;
}

unsigned int Subset::first() const
{
	for (unsigned int i = 0; i < u_size_array; ++i)
	{
		if (ptr_array[i] != 0) {
			uint64_t x = 1;
			unsigned int m = 0;
			while ((ptr_array[i] & x) == 0) {
				x <<= 1;
				++m;
			}
			return (i << 6) + m;
		}
	}
	return ~0u;
}

unsigned int Subset::first_next(unsigned int y) const
{
	++y;
	unsigned int i = (y >> 6);
	if (i >= u_size_array) return ~0u;
	auto x = ptr_array[i];
	unsigned int m = (y & 0x3F);
	x >>= m;
	for (;;) {
		if (x != 0) {
			while ((x & 1) == 0)
			{
				x >>= 1;
				++m;
			}
			return (i << 6) + m;
		}
		else {
			++i;
			if (i != u_size_array) {
				m = 0;
				x = ptr_array[i];
			}
			else return ~0;
		}
	}
}


bool Subset::contains(Subset const & s) const
{
	unsigned int const bound = min(u_size_array, s.u_size_array);
	unsigned int i = 0;
	for (; i < bound; ++i)
	{
		if ((ptr_array[i] & s.ptr_array[i]) != s.ptr_array[i]) return false;
	}
	for (; i < s.u_size_array; ++i)
	{
		if (s.ptr_array[i] != 0) return false;
	}
	return true;
}



void Subset::applyTransform(map<unsigned int, unsigned int> const & my_map)
{
	Subset tmp;
	for (unsigned int i = 0; i < u_size_array; ++i)
	{
		if (ptr_array[i] == 0) continue;
		auto x = ptr_array[i];
		unsigned int n = 64*i;
		do
		{
			while ((x & 1) == 0) {
				x >>= 1;
				++n;
			}
			tmp.emplace(my_map.at(n));
			x >>= 1;
			++n;
		} while (x != 0);
	}
	*this = std::move(tmp);
}

bool disjoint(Subset const & s1, Subset const & s2)
{
	unsigned int const bound = min(s1.u_size_array, s2.u_size_array);
	for (unsigned int i = 0; i < bound; ++i)
	{
		if ((s1.ptr_array[i] & s2.ptr_array[i]) != 0) return false;
	}
	return true;
}
