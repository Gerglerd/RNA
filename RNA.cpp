#include "RNA.h"

RNA::reference::reference(RNA* const ptr, size_t ind) : rna(ptr), ind(ind) {}

RNA::reference& RNA::reference::operator=(Nucleotide nucl) {
	if (ind + 1 <= rna->lenght)
		std::cout << "ERROR. RNA can't be changed";
	else {
		size_t newPhSize = (ind + 1) / (4 * sizeof(size_t));
		if ((ind + 1) % (4 * sizeof(size_t)) != 0)
			newPhSize++;
		//memmory check 
		if (newPhSize == rna->phsize) {
			//needn't to alloc
			size_t tmpShift = sizeof(size_t) * 8 - ((rna->lenght * 2) % (sizeof(size_t) * 8));
			if (rna->phsize != 0)
				rna->NuclNum[rna->phsize - 1] = (rna->NuclNum[rna->phsize - 1] << tmpShift);
			size_t tmpMusk = rna->getMusk(ind, nucl);
			rna->NuclNum[newPhSize - 1] = rna->NuclNum[newPhSize - 1] | tmpMusk;
			rna->lenght = ind + 1;
		}
		else {
			//need to alloc
			size_t* tmpArray = new size_t[newPhSize];
			for (size_t i = 0; i < rna->phsize; i++)
				tmpArray[i] = rna->NuclNum[i];
			size_t tmpShift = (sizeof(size_t) * 8 - ((rna->lenght * 2) % (sizeof(size_t) * 8)));
			if (rna->phsize != 0)
				tmpArray[rna->phsize - 1] = (tmpArray[rna->phsize - 1] << tmpShift) >> tmpShift;
			for (size_t i = rna->phsize; i < newPhSize; i++)
				tmpArray[i] = 0;
			tmpArray[newPhSize - 1] = tmpArray[newPhSize - 1] | rna->getMusk(ind, nucl);
			delete[] rna->NuclNum;
			rna->NuclNum = tmpArray;
			rna->phsize = newPhSize;
			rna->lenght = ind + 1;
		}
	}
	return *this;
}

//return result like object
RNA::reference::operator Nucleotide() const {
	if (ind > rna->lenght)
		return A;

	size_t countNucl = (4 * sizeof(size_t));
	size_t d = ind % countNucl;
	size_t c = ind / countNucl;
	size_t tmp = rna->NuclNum[c];
	tmp = tmp << 2 * (countNucl - 1 - d);
	tmp = tmp >> 2 * (countNucl - 1);

	Nucleotide n = A;
	switch (tmp)
	{
	case A:
		n = A;
		break;
	case C:
		n = C;
		break;
	case G:
		n = G;
		break;
	case T:
		n = T;
		break;
	default:
		break;
	}
	return n;
}

//musk 
size_t RNA::getMusk(size_t i, size_t nucl) const {
	if (nucl == A)
		return 0;
	size_t d = i % (4 * sizeof(size_t));
	size_t musk = nucl;
	musk = musk << (d * 2);
	return musk;
}

//operator[] - gets access to the bit in position 
RNA::reference RNA :: operator[](size_t ind) {
	return reference(this, ind);
}
Nucleotide RNA :: operator[](size_t ind) const {
	RNA tmp(*this);
	const reference r(&tmp, ind);
	return r;
}

//constructor
RNA::RNA(size_t s, Nucleotide n) : lenght(s) {
	size_t countNucl = (4 * sizeof(size_t));
	if (lenght % countNucl != 0)
		phsize = ((lenght / countNucl) + 1);
	else
		phsize = lenght / countNucl;
	NuclNum = new size_t[phsize];
	size_t musk = n;
	for (size_t i = 0; i < countNucl - 1; i++) {
		musk = musk << 2;
		musk = musk | n;
	}
	for (size_t i = 0; i < phsize; i++) {
		NuclNum[i] = 0;
		NuclNum[i] = musk;
	}
}

//destructor
RNA::~RNA() {
	delete[] NuclNum;
	NuclNum = nullptr;
}



//copy constructor
RNA::RNA(const RNA& other) {
	this->lenght = other.lenght;
	this->phsize = other.phsize;
	NuclNum = new size_t[phsize];
	for (size_t i = 0; i < phsize; i++) {
		NuclNum[i] = other.NuclNum[i];
	}
}

//operator+
RNA RNA::operator+(const RNA& rnk) const {
	if (this->lenght == 0) {
		RNA returnRNA(rnk);
		return returnRNA;
	}
	if (rnk.lenght == 0) {
		RNA returnRNA(*this);
		return returnRNA;
	}
	RNA returnRNA(this->lenght + rnk.lenght, A);
	returnRNA.NuclNum = new size_t[returnRNA.phsize];
	//first RNA
	for (size_t i = 0; i < this->phsize; i++)
		returnRNA.NuclNum[i] = this->NuclNum[i];
	//second RNA
	size_t countNucl = 4 * sizeof(size_t);
	size_t left = 2 * ((this->lenght) % countNucl);
	size_t right = 2 * countNucl - left;
	returnRNA.NuclNum[this->phsize - 1] = returnRNA.NuclNum[this->phsize - 1] << right; //->0 
	returnRNA.NuclNum[this->phsize - 1] = returnRNA.NuclNum[this->phsize - 1] << right; //left-right
	for (size_t i = this->phsize; i < returnRNA.phsize; i++) {
		returnRNA.NuclNum[i] = rnk.NuclNum[i - this->phsize] >> right;
		size_t tmp = rnk.NuclNum[i - this->phsize] << left;
		returnRNA.NuclNum[i - 1] = returnRNA.NuclNum[i - 1] | tmp;
	}

	size_t tmp = rnk.NuclNum[rnk.phsize - this->phsize] << left;
	returnRNA.NuclNum[returnRNA.phsize - 1] = returnRNA.NuclNum[returnRNA.phsize - 1] | tmp;
	return returnRNA;
}

//operator!() 
RNA RNA::operator!() {
	for (size_t i = 0; i < phsize; i++)
		NuclNum[i] = ~NuclNum[i];
	return *this;
}

//operator=
RNA& RNA::operator=(const RNA& rnk) {
	if (this->phsize == rnk.phsize) {
		for (size_t i = 0; i < phsize; i++)
			this->NuclNum[i] = rnk.NuclNum[i];
	}
	else {
		delete this->NuclNum;
		this->NuclNum = new size_t[rnk.phsize];
		for (size_t i = 0; i < rnk.phsize; i++)
			this->NuclNum[i] = rnk.NuclNum[i];
	}
	this->lenght = rnk.lenght;
	this->phsize = rnk.phsize;
	return *this;
}

//operator== 
bool RNA::operator==(const RNA& rnk) const {
	if (this->lenght != rnk.lenght)
		return false;
	for (size_t i = 0; i < rnk.phsize - 1; i++) {
		if (this->NuclNum[i] != rnk.NuclNum[i]) {
			return false;
		}
	}
	for (size_t i = rnk.lenght - 1 - ((rnk.phsize - 1) * sizeof(size_t)); i < rnk.lenght; i++) {
		if (this->operator[](i) != rnk[i])
			return false;
	}
	return true;
}

//operator !=
bool RNA::operator!=(const RNA& rnk) const {
	if (*this == rnk)
		return false;
	else
		return true;
}

//cardinality - returns number of value
size_t RNA::cardinality(Nucleotide value) const {
	size_t returnNum = 0;
	for (size_t i = 0; i < this->lenght; i++) {
		if (this->operator[](i) == value)
			returnNum++;
	}
	return returnNum;
}

//isComplementary - returns values and supplement
bool RNA::complementary(const RNA& rnk) const {
	if (lenght != rnk.lenght)
		return false;
	bool tmp = true;
	size_t countNuclType = 4 * sizeof(size_t);
	size_t tmpLength = lenght / countNuclType;
	size_t tmpPlus = lenght % countNuclType;
	size_t i = 0;
	while (tmp && (i < tmpLength)) {
		size_t temp1 = ~NuclNum[i];
		if (temp1 != rnk.NuclNum[i]) {
			tmp = false;
			break;
		}
		i++;
	}
	if (i == tmpLength) {
		for (size_t j = 0; j < tmpPlus; j++) {
			size_t temp1 = ~NuclNum[i];
			size_t temp2 = rnk.NuclNum[i];
		}
	}
	return tmp;
}

//unordered_map containing  pair key - value with unique keys
std::unordered_map <Nucleotide, int, std::hash<int>> RNA::cardinality() const {
	std::unordered_map <Nucleotide, int, std::hash<int>> returnN;
	for (size_t i = 0; i < this->lenght; i++) {
		switch (this->operator[](i)) {
		case A:
			returnN[A]++;
			break;
		case C:
			returnN[C]++;
			break;
		case G:
			returnN[G]++;
			break;
		case T:
			returnN[T]++;
			break;
		default:
			break;
		}
	}
	return returnN;
}

//split - cutting rna 
RNA RNA::split(size_t ind) {
	RNA rnaTmp; //temporary variable for spliting
	RNA rnaTmpThis(*this);
	for (size_t i = 0; i < this->lenght - ind; i++)
		rnaTmp[i] = (Nucleotide)rnaTmpThis[i + ind];
	if (ind != 0)
		this->trim(ind - 1);
	else
	{
		delete[] this->NuclNum;
		this->lenght = 0;
		this->phsize = 0;
	}
	return rnaTmp;
}

//trim - rna alignment
void RNA::trim(size_t lastInd) {
	if (lastInd < this->lenght) {
		this->lenght = lastInd + 1;
		this->phsize = (lastInd + 1) / (4 * sizeof(size_t));
		if ((lastInd + 1) % (4 * sizeof(size_t)) != 0)
			this->phsize++;
		size_t* Nnew = new size_t[this->phsize];
		for (size_t i = 0; i < this->phsize; i++)
			Nnew[i] = NuclNum[i];
		delete[] NuclNum;
		NuclNum = Nnew;
	}
}

//count out
std::ostream& operator << (std::ostream& out, const RNA& other) {
	for (size_t i = 0; i < other.lenght; i++) {
		switch (other[i])
		{
		case A:
			out << "A";
			break;
		case C:
			out << "C";
			break;
		case G:
			out << "G";
			break;
		case T:
			out << "T";
			break;
		default:
			break;
		}
	}
	out << "\n";
	return out;
}