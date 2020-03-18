#include <iostream>
#include <cmath>
using namespace std;
const int N = 6;

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
	/*
	for (int i = 1; i < N; i++) {
		for (int j = 1; j < N; j++) {
			cout << Kf[i][j] << " ";
			if (j == (N - 1)) {
				cout << fi[i] << endl;
			}				
		}
	}*/
	x = gauss(Kf, fi);
	cout << "Gauss method" << endl;
	for (int i = 1; i < N; i++)
		cout << "x" << i << " = " << x[i] << endl;

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
	
	x = Zeidel(Kf,fi);
	
	for (int i = 1; i < N; i++)
		cout << "x" << i << " = " << x[i] << endl;

	return 0;
}
