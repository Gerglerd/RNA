#ifndef RNA_H
#define RNA_H
#include <iostream>
#include <unordered_map>
#include <string>

enum Nucleotide {
	A = 0b000000000,
	C = 0b000000010,
	G = 0b000000001,
	T = 0b000000011
};

class RNA {
private:
	//создаём класс для того, чтобы находить нужный элемент
	friend class reference;
	class reference {
		RNA* rna;
		size_t ind;
	public:
		reference(RNA* const ptr, size_t ind);
		reference& operator = (Nucleotide nucleotide);
		operator Nucleotide() const; //overload operator
	};
	size_t lenght;
	size_t phsize; //physical size
	size_t* NuclNum;
public:
	RNA() : NuclNum(nullptr), lenght(0), phsize(0) {}
	RNA(size_t size, Nucleotide nucl);

	RNA(const RNA& other); //copy-constructor: memmory alloc

	size_t getMusk(size_t i, size_t n) const;

	RNA operator!();
	bool operator==(const RNA& rnk) const;
	bool operator!=(const RNA& rnk) const;
	RNA& operator= (const RNA& rnk);
	RNA operator+(const RNA& rnk) const;

	bool complementary(const RNA& rnk) const;

	void trim(size_t ind);
	RNA split(size_t index);

	size_t cardinality(Nucleotide value) const; //count value for nucleotide
	std::unordered_map <Nucleotide, int, std::hash<int> > cardinality() const; //count value for all types of trits

	reference operator[](size_t ind);
	Nucleotide operator[] (size_t ind) const;
	friend std::ostream& operator << (std::ostream& out, const RNA& other);
	~RNA();
};

#endif
