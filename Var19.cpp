#include <iostream>
#include <cmath>
using namespace std;
const int N = 101;

double* gauss(double kf[][N], double y[]) {
	double* x, max;
	int k, index;
	const double eps = 0.00001;  // точность
	x = new double[N];
	k = 1;
	while (k < N)
	{
		// Поиск строки с максимальным a[i][k]
		max = abs(kf[k][k]);
		index = k;
		for (int i = k + 1; i < N; i++)
		{
			if (abs(kf[i][k]) > max)
			{
				max = abs(kf[i][k]);
				index = i;
			}
		}
		// Перестановка строк
		if (max < eps)
		{
			// нет ненулевых диагональных элементов
			cout << "Решение получить невозможно из-за нулевого столбца ";
			cout << index << " матрицы A" << endl;
			return 0;
		}
		for (int j = 1; j < N; j++)
		{
			double temp = kf[k][j];
			kf[k][j] = kf[index][j];
			kf[index][j] = temp;
		}
		double tempy = y[k];
		y[k] = y[index];
		y[index] = tempy;
		// Нормализация уравнений
		for (int i = k; i < N; i++)
		{
			double temp = kf[i][k];
			if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
			for (int j = 1; j < N; j++)
				kf[i][j] = kf[i][j] / temp;
			y[i] = y[i] / temp;
			if (i == k)  continue; // уравнение не вычитать само из себя
			for (int j = 1; j < N; j++)
				kf[i][j] = kf[i][j] - kf[k][j];
			y[i] = y[i] - y[k];
		}
		k++;
	}
	// обратная подстановка
	for (k = N - 1; k >= 1; k--)
	{
		x[k] = y[k];
		for (int i = 1; i < k; i++)
			y[i] = y[i] - kf[i][k] * x[k];
	}
	return x;
}

double* Zeidel(double kf[][N], double y[]) {
	double eps = 0.00001;
	double* x;
	x = new double[N];
	double max;
	double counter = 0;
	double temp;
	//начальные Х
	for (int i = 0; i < N; i++) {
		x[i] = 0;
	}
	//Итерации
	do {
		for (int i = 1; i < N; i++) {
			temp = x[i];
			x[i] = y[i];
			if (i != (N - 1)) {
				for (int j = i + 1; j < N; j++) {
					x[i] -= x[j] * kf[i][j];
				}
			}
			if (i != 1) {
				for (int j = i - 1; j > 0; j--) {
					x[i] -= x[j] * kf[i][j];
				}
			}
			x[i] = x[i] / kf[i][i];
			//нахождение точности
			if (i == 1) {
				max = abs(x[i] - temp);
			}
			double l = abs(x[i] - temp);
			if (l > max) {
				max = l;
			}
		}

		counter++;
	} while ((max > eps)&&(counter < 100));
	cout << "\nZeidel method" << endl;
	cout << "Iterations = "<< counter << endl;
	return x;
}
//Функция вычисления нормы матрицы
double Norm(double arr[][N]) {
	double norm;
	double sum[N] = {};
	for (int i = 1; i < N; i++) {
		for (int j = 1; j < N; j++) {
			sum[i] += abs(arr[i][j]);
		}
	}

	norm = sum[1];
	for (int i = 2; i < N; i++) {
		if (norm < sum[i]) {
			norm = sum[i];
		}
	}
	return norm;
}
double Norm2(double arr[N]) {
	double norm;
	norm = abs(arr[1]);
	for (int i = 2; i < N; i++) {
		if (norm < abs(arr[i])) {
			norm = abs(arr[i]);
		}
	}
	return norm;
}
//функция нахождения обратной матрицы
void inversion(double A[][N])
{
	double temp;

	double** E = new double* [N];

	for (int i = 0; i < N; i++)
		E[i] = new double[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			E[i][j] = 0.0;

			if (i == j)
				E[i][j] = 1.0;
		}

	for (int k = 1; k < N; k++)
	{
		temp = A[k][k];

		for (int j = 0; j < N; j++)
		{
			A[k][j] /= temp;
			E[k][j] /= temp;
		}

		for (int i = k + 1; i < N; i++)
		{
			temp = A[i][k];

			for (int j = 1; j < N; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	for (int k = N - 1; k > 1; k--)
	{
		for (int i = k - 1; i >= 1; i--)
		{
			temp = A[i][k];

			for (int j = 1; j < N; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	for (int i = 1; i < N; i++)
		for (int j = 1; j < N; j++)
			A[i][j] = E[i][j];

	for (int i = 1; i < N; i++)
		delete[] E[i];

	delete[] E;
}

//Фунция для обусловленности матрицы
double Obusl(double arr[][N]) {
	double Rarr[N][N];//будет обратной матрицей входной матрицы
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			Rarr[i][j] = arr[i][j];
		}
	}
	inversion(Rarr);
	return (Norm(arr)*Norm(Rarr));
}

int main() {
	double Kf[N][N];
	double fi[N];
	double* x;
	x = new double[N];
	x[0] = 0;
	//Зануление всех коэффициентов
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			Kf[i][j] = 0;
		}
	}
	//Присвоение b=10
	for (int i = 1; i < N-1; i++) {
		Kf[i][i] = 10;
	}
	//Присвоение a=1
	for (int i = 2; i < N-1; i++) {
		Kf[i][i-1] = 1;
	}
	//Присвоение c=1
	for (int i = 1; i < N-1; i++) {
		Kf[i][i+1] = 1;
	}
	//Присвоение p=1
	for (int i = 1; i < N; i++) {
		Kf[N-1][i] = 1;
	}
	//Присвоение fi=i
	fi[0] = 0;
	for (int i = 1; i < N; i++) {
		fi[i] = i;
	}
	//Показать матрицу
	/*for (int i = 1; i < N; i++) {
		for (int j = 1; j < N; j++) {
			cout << Kf[i][j] << " ";
			if (j == (N - 1)) {
				cout << fi[i] << endl;
			}				
		}
	}*/
	//cout << "||A|| = " << Norm(Kf)<<endl;
	
	x = gauss(Kf, fi);
	cout << "Gauss method" << endl;
	for (int i = 1; i < N; i++)
		cout << "x" << i << " = " << x[i] << endl;
	cout << "\nobuslovlennost = " << Obusl(Kf)<<endl<<endl;
	//////////////////////////////////////////////////
	//////////////////////////////////////////////////

	//Зануление всех коэффициентов
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			Kf[i][j] = 0;
		}
	}
	//Присвоение b=10
	for (int i = 1; i < N - 1; i++) {
		Kf[i][i] = 10;
	}
	//Присвоение a=1
	for (int i = 2; i < N - 1; i++) {
		Kf[i][i - 1] = 1;
	}
	//Присвоение c=1
	for (int i = 1; i < N - 1; i++) {
		Kf[i][i + 1] = 1;
	}
	//Присвоение p=1
	for (int i = 1; i < N; i++) {
		Kf[N - 1][i] = 1;
	}
	//Присвоение fi=i
	fi[0] = 0;
	for (int i = 1; i < N; i++) {
		fi[i] = i;
	}
	
	//вектор невязки
	double nev[N] = {};
	double sum;
	for (int i = 1; i < N; i++) {
		sum = 0;
		for (int j = 1; j < N; j++) {
			sum += Kf[i][j] * x[j];
		}
		nev[i] = sum - fi[i];
	}
	cout << "Norma neviazki = " << Norm2(nev) << endl;

	//решение методом гаусса-зейделя
	x = Zeidel(Kf,fi);
	
	for (int i = 1; i < N; i++)
		cout << "x" << i << " = " << x[i] << endl;

	for (int i = 1; i < N; i++) {
		sum = 0;
		for (int j = 1; j < N; j++) {
			sum += Kf[i][j] * x[j];
		}
		nev[i] = sum - fi[i];
	}
	cout << "Norma neviazki = " << Norm2(nev) << endl;
	return 0;
}

