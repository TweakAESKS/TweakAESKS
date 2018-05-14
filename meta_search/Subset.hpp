/*
 * Subset.hpp
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


#ifndef DEF_SUBSET
#define DEF_SUBSET

#include <utility>
#include <vector>
#include <iostream>
#include <algorithm>
#include <map>

class Subset
{
	public:
		Subset() : u_size_array(0), ptr_array(nullptr) {};
		Subset(unsigned int x) : u_size_array((x >> 6) + 1), ptr_array(new uint64_t [u_size_array]())
		{
			auto const bound = u_size_array-1;
			ptr_array[bound] = 1;
			ptr_array[bound] <<= x & 0x3F; // x & 0x3F == x % 64
		};

		Subset(Subset const & s) : u_size_array(s.u_size_array), ptr_array(new uint64_t [u_size_array]) {std::copy_n(s.ptr_array, s.u_size_array, ptr_array);};
		Subset(Subset && s) : u_size_array(s.u_size_array), ptr_array(s.ptr_array) {s.u_size_array = 0;	s.ptr_array = nullptr;};

		~Subset() {delete[] ptr_array;};

		Subset & operator=(Subset const & s) {
			if (u_size_array < s.u_size_array)
			{
				delete[] ptr_array;
				ptr_array = new uint64_t [s.u_size_array];
			}
			std::copy_n(s.ptr_array, s.u_size_array, ptr_array);
			u_size_array = s.u_size_array;
			return *this;
		};

		Subset & operator=(Subset && s) {
			u_size_array = s.u_size_array;
			delete[] ptr_array;
			ptr_array = s.ptr_array;
			s.u_size_array = 0;
			s.ptr_array = nullptr;
			return *this;
		};
		Subset & operator=(unsigned int x);

		Subset & operator+=(Subset const &);
		Subset & operator+=(Subset && s) {
			if (u_size_array < s.u_size_array) std::swap(*this,s);
			for (unsigned int i = 0; i < s.u_size_array; ++i) ptr_array[i] |= s.ptr_array[i];
			return *this;
		};
		Subset & operator+=(unsigned int x) {emplace(x); return *this;};
		Subset & operator-=(Subset const & s) {
			auto const bound = std::min(u_size_array, s.u_size_array);
			for (unsigned int i = 0; i < bound; ++i) ptr_array[i] &= ~s.ptr_array[i];
			return *this;
		};

		Subset & operator-=(unsigned int x) {erase(x); return *this;};
		Subset & operator&=(Subset const & s) {
			if (u_size_array > s.u_size_array) u_size_array = s.u_size_array;
			for (unsigned int i = 0; i < u_size_array; ++i) ptr_array[i] &= s.ptr_array[i];
			return *this;
		};

		Subset & operator^=(Subset const &);
		Subset & operator^=(Subset &&);

		unsigned int size() const;
		bool empty() const;
		void clear() {u_size_array = 0;};

		void emplace(unsigned int x);
		void erase(unsigned int x) {auto const bound = (x >> 6); if (u_size_array > bound) ptr_array[bound] &= ~(uint64_t (1) << (x & 0x3F));};

		bool contains(Subset const &) const;
		bool contains(unsigned int x) const {
			unsigned int const bound = (x >> 6);
			return ((bound < u_size_array) && (((ptr_array[bound] >> (x & 0x3F)) & 1) == 1));
		};

		unsigned int count(unsigned int x) const {return contains(x) ? 1 : 0;};
		bool shareElements(Subset const & s) const {
			unsigned int const bound = std::min(u_size_array, s.u_size_array);
			for (unsigned int i = 0; i < bound; ++i)
			{
				if ((ptr_array[i] & s.ptr_array[i]) != 0) return true;
			}
			return false;
		};
		bool shareElements(unsigned int x) const {return contains(x);};

		std::vector<unsigned int> getElements() const;
	  unsigned int first() const;
		unsigned int first_next(unsigned int) const;

		template <typename T>
		unsigned int applyThenSum(T f) const;

		template <typename T>
		void apply(T f) const;

		void applyTransform(std::map<unsigned int, unsigned int> const &);

		template <typename T>
		void applyThenRemoveIfFalse(T f);

		friend Subset operator+(Subset const & s1, Subset const & s2) {return (s1.u_size_array < s2.u_size_array) ? s1.getUnionWith(s2):s2.getUnionWith(s1);};
		friend Subset operator+(Subset && s1, Subset const & s2) {return (s1.u_size_array < s2.u_size_array) ? s1.getUnionWith(s2):s2.getUnionWith(std::move(s1));};
		friend Subset operator+(Subset const & s1, Subset && s2) {return (s1.u_size_array < s2.u_size_array) ? s1.getUnionWith(std::move(s2)):s2.getUnionWith(s1);};
		friend Subset operator+(Subset && s1, Subset && s2) {return (s1.u_size_array < s2.u_size_array) ? s1.getUnionWith(std::move(s2)):s2.getUnionWith(std::move(s1));};

		friend Subset operator-(Subset const & s1, Subset const & s2) {auto s = s1; s -= s2; return s;};
		friend Subset operator-(Subset && s1, Subset const & s2) {auto s = std::move(s1); s -= s2; return s;};

		friend Subset operator&(Subset const & s1, Subset const & s2) {
			const unsigned int size_ptr = std::min(s1.u_size_array, s2.u_size_array);
			uint64_t * __restrict__ ptr = new uint64_t [size_ptr];
			for (unsigned int i = 0; i < size_ptr; ++i) ptr[i] = s1.ptr_array[i] & s2.ptr_array[i];
			return Subset(size_ptr, ptr);
		};
		friend Subset operator&(Subset && s1, Subset const & s2) {auto s = std::move(s1); s &= s2; return s;};
		friend Subset operator&(Subset const & s1, Subset && s2) {auto s = std::move(s2); s &= s1; return s;};
		friend Subset operator&(Subset && s1, Subset && s2) {auto s = std::move(s1); s &= s2; return s;};

		friend Subset operator^(Subset const & s1, Subset const & s2) {return (s1.u_size_array < s2.u_size_array) ? s1.getXorWith(s2):s2.getXorWith(s1);};
		friend Subset operator^(Subset && s1, Subset const & s2) {return (s1.u_size_array < s2.u_size_array) ? s1.getXorWith(s2):s2.getXorWith(std::move(s1));};
		friend Subset operator^(Subset const & s1, Subset && s2) {return (s1.u_size_array < s2.u_size_array) ? s1.getXorWith(std::move(s2)):s2.getXorWith(s1);};
		friend Subset operator^(Subset && s1, Subset && s2) {return (s1.u_size_array < s2.u_size_array) ? s1.getXorWith(std::move(s2)):s2.getXorWith(std::move(s1));};

		friend void swap(Subset & s1, Subset & s2) {std::swap(s1.u_size_array, s2.u_size_array); std::swap(s1.ptr_array, s2.ptr_array);};

		friend bool disjoint(Subset const & s1, Subset const & s2);

		friend std::ostream & operator<<( std::ostream &flux, Subset const & s) {
			for (unsigned int i = 0; i < s.u_size_array; ++i) {
				auto x = s.ptr_array[i];
				for (unsigned j = 0; j < 64; ++j) flux << ((x >> j) & 1);
			}
			return flux;
		}

		// Relational Operators (for fast comparison)
		friend bool operator<(Subset const & s1, Subset const & s2) {
			if (s1.u_size_array != s2.u_size_array) return (s1.u_size_array < s2.u_size_array);
			for (unsigned int i = 0; i < s1.u_size_array; ++i) {
				if (s1.ptr_array[i] != s2.ptr_array[i]) return (s1.ptr_array[i] < s2.ptr_array[i]);
			}
			return false;
		};
		friend bool operator<=(Subset const & s1, Subset const & s2) {
			if (s1.u_size_array != s2.u_size_array) return (s1.u_size_array < s2.u_size_array);
			for (unsigned int i = 0; i < s1.u_size_array; ++i) {
				if (s1.ptr_array[i] != s2.ptr_array[i]) return (s1.ptr_array[i] < s2.ptr_array[i]);
			}
			return true;
		};
		friend bool operator>(Subset const & s1, Subset const & s2) {return s2 < s1;};
		friend bool operator>=(Subset const & s1, Subset const & s2) {return s2 <= s1;};

		friend bool operator==(Subset const & s1, Subset const & s2) {
			if (s1.u_size_array != s2.u_size_array) return false;
			for (unsigned int i = 0; i < s1.u_size_array; ++i) {
				if (s1.ptr_array[i] != s2.ptr_array[i]) return false;
			}
			return true;
		};
		friend bool operator!=(Subset const & s1, Subset const & s2) {return !(s1 == s2);};

	private:
		unsigned int u_size_array;
		uint64_t * ptr_array;


		Subset(unsigned int size_a, uint64_t * a) : u_size_array(size_a), ptr_array(a) {};

		Subset getUnionWith(Subset const & s) const {
			uint64_t * __restrict__ ptr = new uint64_t [s.u_size_array];
			unsigned int i = 0;
			for (; i < u_size_array; ++i) ptr[i] = ptr_array[i] | s.ptr_array[i];
			for (; i < s.u_size_array; ++i) ptr[i] = s.ptr_array[i];
			return Subset(s.u_size_array, ptr);
		};
		Subset getUnionWith(Subset && s) const {
			for (unsigned int i = 0; i < u_size_array; ++i) s.ptr_array[i] |= ptr_array[i];
			return s;
		};

		Subset getXorWith(Subset const & s) const {
			uint64_t * __restrict__ ptr = new uint64_t [s.u_size_array];
			unsigned int i = 0;
			for (; i < u_size_array; ++i) ptr[i] = s.ptr_array[i] ^ ptr_array[i];
			for (; i < s.u_size_array; ++i) ptr[i] = s.ptr_array[i];
			return Subset(s.u_size_array, ptr);
		};
		Subset getXorWith(Subset && s) const {
			for (unsigned int i = 0; i < u_size_array; ++i) s.ptr_array[i] ^= ptr_array[i];
			return s;
		};

};


template <typename T>
unsigned int Subset::applyThenSum(T f) const
{
	unsigned int res = 0;
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
			res += f(n);
			x >>= 1;
			++n;
		} while (x != 0);
	}
	return res;
}

template <typename T>
void Subset::apply(T f) const
{
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
			f(n);
			x >>= 1;
			++n;
		} while (x != 0);
	}
}

template <typename T>
void Subset::applyThenRemoveIfFalse(T f)
{
	for (unsigned int i = 0; i < u_size_array; ++i)
	{
		if (ptr_array[i] == 0) continue;
		auto x = ptr_array[i];
		unsigned int const n = 64*i;
		unsigned int m = 0;
		do
		{
			while ((x & 1) == 0) {
				x >>= 1;
				++m;
			}
			if (f(n+m) == false) ptr_array[i] ^= (uint64_t (1)) << m;
			x >>= 1;
			++m;
		} while (x != 0);
	}
}


#endif
