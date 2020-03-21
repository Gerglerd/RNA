#include <iostream>
#include "RNA.h"
#include <ctime>

int main() {
	setlocale(LC_ALL, "");
	RNA rna1;
	size_t size = 25;
	RNA rna2(size, G);
	RNA rna3(rna2);
	std::cout  << "rna2: " << rna2 << "\n" << "rna3: " << rna3 << std::endl;
	int i = 0;
	rna1[0] = T;
	rna1[1] = G;
	rna1[2] = C;
	rna1[3] = A;
	rna1[4] = T;
	rna1[5] = T;
	rna1[6] = G;
	rna1[7] = A;
	rna1[8] = A;
	rna1[9] = T;
	rna1[10] = C;
	rna1[11] = G;
	rna1[12] = G;
	rna1[13] = G;
	rna1[14] = T;
	rna1[15] = T;
	rna1[16] = A;
	rna1[17] = T;
	rna1[18] = C;
	rna1[19] = C;
	rna1[20] = A;
	rna1[21] = C;
	rna1[22] = T;
	rna1[23] = C;
	rna1[24] = G;
	rna1[25] = T;
	rna1[26] = T;
	RNA r;
	unsigned int start_time = clock();
	for (size_t i = 0; i < 1000000; i++)
	{
		r[i] = A;
		if (i == 999999)
			std::cout << "the last element: " <<r[i] << std::endl;
	}
	unsigned int end_time = clock();
	unsigned int search_time = end_time - start_time;
	std::cout << "testing time: " << search_time << std::endl;
	std::cout << "rna1: " << rna1 << std::endl;
	std::cout << "rna1 + rna2: " << rna1 + rna2 << std::endl;
	rna2 = rna1.split(5);
	std::cout << "rna2 = rna1.split to 5: " << rna2;
	system("pause");
	return 0;
}