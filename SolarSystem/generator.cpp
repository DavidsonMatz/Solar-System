#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
using namespace std;

vector<float> splitfloat(const string& s, char d, vector<float> vetor) {
	stringstream ss(s);
	string texto;
	while (getline(ss, texto, d)) {
		vetor.push_back(stof(texto));
	}
	return vetor;
}

vector<float> cross(vector<float> a, vector<float> b, vector<float> res) {

	res[0] = a[1] * b[2] - a[2] * b[1];
	res[1] = a[2] * b[0] - a[0] * b[2];
	res[2] = a[0] * b[1] - a[1] * b[0];
	return res;
}

vector<float> normalize(vector<float> a) {

	float l = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
	a[0] = a[0] / l;
	a[1] = a[1] / l;
	a[2] = a[2] / l;
	return a;
}

void transpose(float m[][4], float res[][4]){
	for (int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
			res[i][j] = m[j][i];
		}
	}
}

void plano(float comprimento, float largura, const char* file) {
	comprimento /= 2;
	largura /= 2;
	ofstream f;
	f.open(file);

	f
		<< comprimento << "," << 0.0f << "," << largura << "\n"
		<< 0.0f << "," << 1.0f << "," << 0.0f << "\n"
		<< 1.0f << "," << 1.0f << "\n"
		<< -comprimento << "," << 0.0f << "," << -largura << "\n"
		<< 0.0f << "," << 1.0f << "," << 0.0f << "\n"
		<< 0.0f << "," << 0.0f << "\n"
		<< -comprimento << "," << 0.0f << "," << largura << "\n"
		<< 0.0f << "," << 1.0f << "," << 0.0f << "\n"
		<< 0.0f << "," << 1.0f << "\n"

		<< comprimento << "," << 0.0f << "," << largura << "\n"
		<< 0.0f << "," << 1.0f << "," << 0.0f << "\n"
		<< 1.0f << "," << 1.0f << "\n"
		<< comprimento << "," << 0.0f << "," << -largura << "\n"
		<< 0.0f << "," << 1.0f << "," << 0.0f << "\n"
		<< 1.0f << "," << 0.0f << "\n"
		<< -comprimento << "," << 0.0f << "," << -largura << "\n"
		<< 0.0f << "," << 1.0f << "," << 0.0f << "\n"
		<< 0.0f << "," << 0.0f << "\n"

		<< comprimento << "," << 0.0f << "," << largura << "\n"
		<< 0.0f << "," << -1.0f << "," << 0.0f << "\n"
		<< 1.0f << "," << 1.0f << "\n"
		<< -comprimento << "," << 0.0f << "," << largura << "\n"
		<< 0.0f << "," << -1.0f << "," << 0.0f << "\n"
		<< 0.0f << "," << 1.0f << "\n"
		<< -comprimento << "," << 0.0f << "," << -largura << "\n"
		<< 0.0f << "," << -1.0f << "," << 0.0f << "\n"
		<< 0.0f << "," << 0.0f << "\n"


		<< comprimento << "," << 0.0f << "," << largura << "\n"
		<< 0.0f << "," << -1.0f << "," << 0.0f << "\n"
		<< 1.0f << "," << 1.0f << "\n"
		<< -comprimento << "," << 0.0f << "," << -largura << "\n"
		<< 0.0f << "," << -1.0f << "," << 0.0f << "\n"
		<< 0.0f << "," << 0.0f << "\n"
		<< comprimento << "," << 0.0f << "," << -largura << "\n"
		<< 0.0f << "," << -1.0f << "," << 0.0f << "\n"
		<< 1.0f << "," << 0.0f << "\n";

	f.close();
}

void cilindro(float raio, float altura, float slices, const char* file) {
	int i;
	float px1 = 0.0f;
	float pz1 = raio;
	float px2, pz2;
	float alpha;
	ofstream f;
	f.open(file);

	for (int i = 1; i <= slices; i++) {
		alpha = (2.00f * M_PI * float(i)) / float(slices);
		px2 = raio * sin(alpha);
		pz2 = raio * cos(alpha);

		// BASE DE BAIXO
		f << 0.0f << "," << 0.0f << "," << 0.0f << "\n"
			<< 0.0f << "," << -1.0f << "," << 0.0f << "\n"
			<< 0.5f << "," << 0.5f << "\n"
			<< px2 << "," << 0.0f << "," << pz2 << "\n"
			<< 0.0f << "," << -1.0f << "," << 0.0f << "\n"
			<< (px2 / raio) * 0.5f + 0.5f << "," << (pz2 / raio) * 0.5f + 0.5f << "\n"
			<< px1 << "," << 0.0f << "," << pz1 << "\n"
			<< 0.0f << "," << -1.0f << "," << 0.0f << "\n"
			<< (px1 / raio) * 0.5f + 0.5f << "," << (pz1 / raio) * 0.5f + 0.5f << "\n"

			// BASE DE CIMA
			<< 0.0f << "," << altura << "," << 0.0f << "\n"
			<< 0.0f << "," << 1.0f << "," << 0.0f << "\n"
			<< 0.5f << "," <<  0.5f << "\n"
			<< px1 << "," << altura << "," << pz1 << "\n"
			<< 0.0f << "," << 1.0f << "," << 0.0f << "\n"
			<< (px1 / raio) * 0.5f + 0.5f << "," << (pz1 / raio) * 0.5f + 0.5f << "\n"
			<< px2 << "," << altura << "," << pz2 << "\n"
			<< 0.0f << "," << 1.0f << "," << 0.0f << "\n"
			<< (px2 / raio) * 0.5f + 0.5f << "," << (pz2 / raio) * 0.5f + 0.5f << "\n"

			//LADOS
			<< px1 << "," << 0.0f << "," << pz1 << "\n"
			<< px1 / raio << "," << 0.0f << "," << pz1 / raio << "\n" 
			<< (i - 1) / float(slices) << "," << 0.0f << "\n"
			<< px2 << "," << altura << "," << pz2 << "\n"
			<< px2 / raio << "," << 0.0f << "," << pz2 / raio << "\n"
			<< i / float(slices) << "," << 1.0f << "\n"
			<< px1 << "," << altura << "," << pz1 << "\n"
			<< px1 / raio << "," << 0.0f << "," << pz1 / raio << "\n"
			<< (i - 1) / float(slices) << "," << 1.0f << "\n"

			<< px1 << "," << 0.0f << "," << pz1 << "\n"
			<< px1 / raio << "," << 0.0f << "," << pz1 / raio << "\n"
			<< (i - 1) / float(slices) << "," << 0.0f << "\n"
			<< px2 << "," << 0.0f << "," << pz2 << "\n"
			<< px2 / raio << "," << 0.0f << "," << pz2 / raio << "\n"
			<< i / float(slices) << "," << 0.0f << "\n"
			<< px2 << "," << altura << "," << pz2 << "\n"
			<< px2 / raio << "," << 0.0f << "," << pz2 / raio << "\n"
			<< i / float(slices) << "," << 1.0f << "\n";


		px1 = px2;
		pz1 = pz2;
	}
	f.close();
}

void torus(float raio_int, float raio_ext, float slices, float stacks, const char* file) {
    ofstream f;
	f.open(file);
    float alpha = 0, beta = 0;
    float incr_a = (2*M_PI)/slices, incr_b = (2*M_PI)/stacks;
    float px1, py1, pz1, px2, py2, pz2, px3, py3, pz3, px4, py4, pz4;
	float px1n, py1n, pz1n, px2n, py2n, pz2n, px3n, py3n, pz3n, px4n, py4n, pz4n;
	float t1, t2;
	int half_stack = stacks / 2;
    for (int i = 0; i < slices; i++) {
        for(int j = 0; j < stacks; j++) {

            px1 = cos(alpha) * (raio_ext + raio_int * cos(beta));
            py1 = sin(alpha) * (raio_ext + raio_int * cos(beta));
            pz1 = raio_int * sin(beta);

            px2 = cos(alpha + incr_a) * (raio_ext + raio_int * cos(beta));
            py2 = sin(alpha + incr_a) * (raio_ext + raio_int * cos(beta));
            pz2 = raio_int * sin(beta);

            px3 = cos(alpha + incr_a) * (raio_ext + raio_int * cos(beta + incr_b));
            py3 = sin(alpha + incr_a) * (raio_ext + raio_int * cos(beta + incr_b));
            pz3 = raio_int * sin(beta + incr_b);

            px4 =  cos(alpha) * (raio_ext + raio_int * cos(beta + incr_b));
            py4 =  sin(alpha) * (raio_ext + raio_int * cos(beta + incr_b));
            pz4 =  raio_int * sin(beta + incr_b);


			px1n = cos(alpha) * cos(beta);
			py1n = sin(alpha) * cos(beta);
			pz1n = sin(beta);

			px2n = cos(alpha + incr_a) * cos(beta);
			py2n = sin(alpha + incr_a) * cos(beta);
			pz2n = sin(beta);

			px3n = cos(alpha + incr_a) * cos(beta + incr_b);
			py3n = sin(alpha + incr_a) * cos(beta + incr_b);
			pz3n = sin(beta + incr_b);

			px4n = cos(alpha) * cos(beta + incr_b);
			py4n = sin(alpha) * cos(beta + incr_b);
			pz4n = sin(beta + incr_b);
			
			if (j > half_stack) {
				t1 = 1 - ((j / (stacks / 2)) - 1);
				t2 = 1 - (((j + 1) / (stacks / 2)) - 1);
			}
			else {
				t1 = j / (stacks / 2);
				t2 = (j + 1) / (stacks / 2);
			}

			f
				<< px1 << "," << py1 << "," << pz1 << "\n"
				<< px1n << "," << py1n << "," << pz1n << "\n"
				<< 0.0f << "," << t1 << "\n"
				<< px2 << "," << py2 << "," << pz2 << "\n"
				<< px2n << "," << py2n << "," << pz2n << "\n"
				<< 1.0f << "," << t1 << "\n"
				<< px3 << "," << py3 << "," << pz3 << "\n"
				<< px3n << "," << py3n << "," << pz3n << "\n"
				<< 1.0f << "," << t2 << "\n"

				<< px3 << "," << py3 << "," << pz3 << "\n"
				<< px3n << "," << py3n << "," << pz3n << "\n"
				<< 1.0f << "," << t2 << "\n"
				<< px4 << "," << py4 << "," << pz4 << "\n"
				<< px4n << "," << py4n << "," << pz4n << "\n"
				<< 0.0f << "," << t2 << "\n"
				<< px1 << "," << py1 << "," << pz1 << "\n"
				<< px1n << "," << py1n << "," << pz1n << "\n"
				<< 0.0f << "," << t1 << "\n";

			beta = incr_b * (j + 1);
        }
        alpha = incr_a * (i + 1);
    }
    f.close();
}

void esfera(float raio, float slices, float stacks, const char* file) {
	float alpha, beta;
	float px2, py2, pz2, px3, py3, pz3, px4, pz4, pxdecima2, pzdecima2;
	ofstream f;
	f.open(file);
	float px1 = 0.0f;
	float py1 = -raio;
	float pz1 = 0.0f;
	for (int j = 0; j < stacks; j++) {
		beta = (-M_PI / 2 + M_PI * float(j) / float(stacks));
		for (int i = 0; i < slices; i++) {
			alpha = (2.0f * M_PI * float(i) / float(slices));
			px1 = raio * cos(beta) * sin(alpha);
			py1 = raio * sin(beta);
			pz1 = raio * cos(beta) * cos(alpha);
			px2 = raio * cos(beta) * sin(alpha + 2 * M_PI / slices);
			pz2 = raio * cos(beta) * cos(alpha + 2 * M_PI / slices);
			px3 = raio * cos(beta + M_PI / stacks) * sin(alpha + 2 * M_PI / slices);
			py3 = raio * sin(beta + M_PI / stacks);
			pz3 = raio * cos(beta + M_PI / stacks) * cos(alpha + 2 * M_PI / slices);
			px4 = raio * cos(beta + M_PI / stacks) * sin(alpha);
			pz4 = raio * cos(beta + M_PI / stacks) * cos(alpha);


			f
				<< px1 << "," << py1 << "," << pz1 << "\n"
				<< px1 / raio << "," << py1 / raio << "," << pz1 / raio << "\n"
				<< float(i/slices) << "," << float(j/stacks) << "\n"
				<< px2 << "," << py1 << "," << pz2 << "\n"
				<< px2 / raio << "," << py1 / raio << "," << pz2 / raio << "\n"
				<< float((i+1) / slices) << "," << float(j / stacks) << "\n"
				<< px3 << "," << py3 << "," << pz3 << "\n"
				<< px3 / raio << "," << py3 / raio << "," << pz3 / raio << "\n"
				<< float((i + 1) / slices) << "," << float((j+1) / stacks) << "\n"

				<< px1 << "," << py1 << "," << pz1 << "\n"
				<< px1 / raio << "," << py1 / raio << "," << pz1 / raio << "\n"
				<< float(i / slices) << "," << float(j / stacks) << "\n"
				<< px3 << "," << py3 << "," << pz3 << "\n"
				<< px3 / raio << "," << py3 / raio << "," << pz3 / raio << "\n"
				<< float((i + 1) / slices) << "," << float((j + 1)) / stacks << "\n"
				<< px4 << "," << py3 << "," << pz4 << "\n"
				<< px4 / raio << "," << py3 / raio << "," << pz4 / raio << "\n"
				<< float(i / slices) << "," << float((j + 1) / stacks) << "\n"
				;
		}
	}
	f.close();
}

void cone(float raio, float altura, float slices, float stacks, const char* file) {
	float px2, pz2, pxdecima2, pzdecima2;
	float alpha, alpha2;
	ofstream f;
	float angle = atan(raio / altura);
	f.open(file);
	for (int j = stacks; j >= 1; j--) {
		float px1 = 0.0f;
		float pz1 = raio * j / stacks;
		float pxdecima1 = 0.0f;
		float pzdecima1 = (raio / stacks) * (j - 1);
		for (int i = 1; i <= slices; i++) {
			alpha = (2.00f * M_PI * float(i) / float(slices));
			alpha2 = (2.00f * M_PI * float(i-1) / float(slices));
			px2 = raio * sin(alpha) * j / stacks;
			pz2 = raio * cos(alpha) * j / stacks;
			pxdecima2 = (raio / stacks) * (j - 1) * sin(alpha);
			pzdecima2 = (raio / stacks) * (j - 1) * cos(alpha);

			// BASE DE BAIXO
			if (j == stacks) {
				f
					<< 0.0f << "," << 0.0f << "," << 0.0f << "\n"
					<< 0.0f << "," << -1.0f << "," << 0.0f << "\n"
					<< 0.5f << ","  << 0.5f << "\n"
					<< px2 << "," << 0.0f << "," << pz2 << "\n"
					<< 0.0f << "," << -1.0f << "," << 0.0f << "\n"
					<< (px2 / raio) * 0.5f + 0.5f << "," << (pz2 / raio) * 0.5f + 0.5f << "\n"
					<< px1 << "," << 0.0f << "," << pz1 << "\n"
					<< 0.0f << "," << -1.0f << "," << 0.0f << "\n"
					<< (px1 / raio) * 0.5f + 0.5f << "," << (pz1 / raio) * 0.5f + 0.5f << "\n";
			}
			//LADOS
			if (j == 1) {
				f
					<< px1 << "," << (stacks - j) * (altura / stacks) << "," << pz1 << "\n"
					<< cos(angle) * sin(alpha2) << "," << sin(angle) << "," << cos(angle) * cos(alpha2) << "\n"
					<< 0.0f << "," << 0.0f << "\n"
					<< px2 << "," << (stacks - j) * (altura / stacks) << "," << pz2 << "\n"
					<< cos(angle) * sin(alpha) << "," << sin(angle) << "," << cos(angle) * cos(alpha) << "\n"
					<< 1.0f << "," << 0.0f << "\n"
					<< pxdecima2 << "," << (stacks - j + 1) * (altura / stacks) << "," << pzdecima2 << "\n"
					<< 0.0f << "," << 1.0f << "," << 0.0f << "\n"
					<< 1.0f << "," << 1.0f << "\n"

					<< px1 << "," << (stacks - j) * (altura / stacks) << "," << pz1 << "\n"
					<< cos(angle) * sin(alpha2) << "," << sin(angle) << "," << cos(angle) * cos(alpha2) << "\n"
					<< 0.0f << "," << 0.0f << "\n"
					<< pxdecima2 << "," << (stacks - j + 1) * (altura / stacks) << "," << pzdecima2 << "\n"
					<< 0.0f << "," << 1.0f << "," << 0.0f << "\n"
					<< 1.0f << "," << 1.0f << "\n"
					<< pxdecima1 << "," << (stacks - j + 1) * (altura / stacks) << "," << pzdecima1 << "\n"
					<< 0.0f << "," << 1.0f << "," << 0.0f << "\n"
					<< 0.0f << "," << 1.0f << "\n";
			}
			else {
				f
					<< px1 << "," << (stacks - j) * (altura / stacks) << "," << pz1 << "\n"
					<< cos(angle) * sin(alpha2) << "," << sin(angle) << "," << cos(angle) * cos(alpha2) << "\n"
					<< 0.0f << "," << 0.0f << "\n"
					<< px2 << "," << (stacks - j) * (altura / stacks) << "," << pz2 << "\n"
					<< cos(angle) * sin(alpha) << "," << sin(angle) << "," << cos(angle) * cos(alpha) << "\n"
					<< 1.0f << "," << 0.0f << "\n"
					<< pxdecima2 << "," << (stacks - j + 1) * (altura / stacks) << "," << pzdecima2 << "\n"
					<< cos(angle) * sin(alpha) << "," << sin(angle) << "," << cos(angle) * cos(alpha) << "\n"
					<< 1.0f << "," << 1.0f << "\n"

					<< px1 << "," << (stacks - j) * (altura / stacks) << "," << pz1 << "\n"
					<< cos(angle) * sin(alpha2) << "," << sin(angle) << "," << cos(angle) * cos(alpha2) << "\n"
					<< 0.0f << "," << 0.0f << "\n"
					<< pxdecima2 << "," << (stacks - j + 1) * (altura / stacks) << "," << pzdecima2 << "\n"
					<< cos(angle) * sin(alpha) << "," << sin(angle) << "," << cos(angle) * cos(alpha) << "\n"
					<< 1.0f << "," << 1.0f << "\n"
					<< pxdecima1 << "," << (stacks - j + 1) * (altura / stacks) << "," << pzdecima1 << "\n"
					<< cos(angle) * sin(alpha2) << "," << sin(angle) << "," << cos(angle) * cos(alpha2) << "\n"
					<< 0.0f << "," << 1.0f << "\n";
			}
				/*
			f
				<< px1 << "," << (stacks - j) * (altura / stacks) << "," << pz1 << "\n"
				<< cos(angle) * px1 / raio << "," << sin(angle) << "," << cos(angle) * pz1/ raio << "\n"
				<< 0.0f << "," << 0.0f << "\n"
				<< px2 << "," << (stacks - j) * (altura / stacks) << "," << pz2 << "\n"
				<< cos(angle) * px2 / raio << "," << sin(angle) << "," << cos(angle) * pz2 / raio << "\n"
				<< 1.0f << "," << 0.0f << "\n"
				<< pxdecima2 << "," << (stacks - j + 1) * (altura / stacks) << "," << pzdecima2 << "\n"
				<< cos(angle) * pxdecima2 / raio << "," << sin(angle) << "," << cos(angle) * pzdecima2 / raio << "\n"
				<< 1.0f << "," << 1.0f << "\n"

				<< px1 << "," << (stacks - j) * (altura / stacks) << "," << pz1 << "\n"
				<< cos(angle) * px1 / raio << "," << sin(angle) << "," << cos(angle) * pz1 / raio << "\n"
				<< 0.0f << "," << 0.0f << "\n"
				<< pxdecima2 << "," << (stacks - j + 1) * (altura / stacks) << "," << pzdecima2 << "\n"
				<< cos(angle) * pxdecima2 / raio << "," << sin(angle) << "," << cos(angle) * pzdecima2 / raio << "\n"
				<< 1.0f << "," << 1.0f << "\n"
				<< pxdecima1 << "," << (stacks - j + 1) * (altura / stacks) << "," << pzdecima1 << "\n"
				<< cos(angle) * pxdecima1 / raio << "," << sin(angle) << "," << cos(angle) * pzdecima1 / raio << "\n"
				<< 0.0f << "," << 1.0f << "\n";
				*/

			pxdecima1 = pxdecima2;
			pzdecima1 = pzdecima2;
			px1 = px2;
			pz1 = pz2;
		}

	}
	f.close();
}

void caixa(float c, float l, float a, int div, const char* file) {
	c /= 2;
	l /= 2;
	float x = a / 2;
	ofstream f;

	f.open(file);
	for (int i = 0; i < div + 1; i++) {
		for (int j = 0; j < div + 1; j++) {
			f
				//base de baixo
				<< -c + (2 * c * i) / float(div + 1) << "," << 0.0f - x << "," << -l + (2 * l * j) / float(div + 1) << "\n"
				<< 0.0f << "," << -1.0f << "," << 0.0f << "\n"
				<< 0.0f << "," << 0.0f << "\n"
				<< -c + (2 * c * (i + 1)) / float(div + 1) << "," << 0.0f - x << "," << -l + (2 * l * (j + 1)) / float(div + 1) << "\n"
				<< 0.0f << "," << -1.0f << "," << 0.0f << "\n"
				<< 1.0f << "," << 1.0f << "\n"
				<< -c + (2 * c * i) / float(div + 1) << "," << 0.0f - x << "," << -l + (2 * l * (j + 1)) / float(div + 1) << "\n"
				<< 0.0f << "," << -1.0f << "," << 0.0f << "\n"
				<< 0.0f << "," << 1.0f << "\n"

				<< -c + (2 * c * i) / float(div + 1) << "," << 0.0f - x << "," << -l + (2 * l * j) / float(div + 1) << "\n"
				<< 0.0f << "," << -1.0f << "," << 0.0f << "\n"
				<< 0.0f << "," << 0.0f << "\n"
				<< -c + (2 * c * (i + 1)) / float(div + 1) << "," << 0.0f - x << "," << -l + (2 * l * j) / float(div + 1) << "\n"
				<< 0.0f << "," << -1.0f << "," << 0.0f << "\n"
				<< 1.0f << "," << 0.0f << "\n"
				<< -c + (2 * c * (i + 1)) / float(div + 1) << "," << 0.0f - x << "," << -l + (2 * l * (j + 1)) / float(div + 1) << "\n"
				<< 0.0f << "," << -1.0f << "," << 0.0f << "\n"
				<< 1.0f << "," << 1.0f << "\n"

				//base de cima
				<< -c + (2 * c * i) / float(div + 1) << "," << a - x << "," << -l + (2 * l * j) / float(div + 1) << "\n"
				<< 0.0f << "," << 1.0f << "," << 0.0f << "\n"
				<< 0.0f << "," << 1.0f << "\n"
				<< -c + (2 * c * i) / float(div + 1) << "," << a - x << "," << -l + (2 * l * (j + 1)) / float(div + 1) << "\n"
				<< 0.0f << "," << 1.0f << "," << 0.0f << "\n"
				<< 0.0f << "," << 0.0f << "\n"
				<< -c + (2 * c * (i + 1)) / float(div + 1) << "," << a - x << "," << -l + (2 * l * (j + 1)) / float(div + 1) << "\n"
				<< 0.0f << "," << 1.0f << "," << 0.0f << "\n"
				<< 1.0f << "," << 0.0f << "\n"

				<< -c + (2 * c * i) / float(div + 1) << "," << a - x << "," << -l + (2 * l * j) / float(div + 1) << "\n"
				<< 0.0f << "," << 1.0f << "," << 0.0f << "\n"
				<< 0.0f << "," << 1.0f << "\n"
				<< -c + (2 * c * (i + 1)) / float(div + 1) << "," << a - x << "," << -l + (2 * l * (j + 1)) / float(div + 1) << "\n"
				<< 0.0f << "," << 1.0f << "," << 0.0f << "\n"
				<< 1.0f << "," << 0.0f << "\n"
				<< -c + (2 * c * (i + 1)) / float(div + 1) << "," << a - x << "," << -l + (2 * l * j) / float(div + 1) << "\n"
				<< 0.0f << "," << 1.0f << "," << 0.0f << "\n"
				<< 1.0f << "," << 1.0f << "\n"

				//lado de frente( L sempre igual)
				<< -c + (2 * c * i) / float(div + 1) << "," << a * j / float(div + 1) - x << "," << l << "\n"
				<< 0.0f << "," << 0.0f << "," << 1.0f << "\n"
				<< 0.0f << "," << 0.0f << "\n"
				<< -c + (2 * c * (i + 1)) / float(div + 1) << "," << a * j / float(div + 1) - x << "," << l << "\n"
				<< 0.0f << "," << 0.0f << "," << 1.0f << "\n"
				<< 1.0f << "," << 0.0f << "\n"
				<< -c + (2 * c * (i + 1)) / float(div + 1) << "," << a * (j + 1) / float(div + 1) - x << "," << l << "\n"
				<< 0.0f << "," << 0.0f << "," << 1.0f << "\n"
				<< 1.0f << "," << 1.0f << "\n"

				<< -c + (2 * c * i) / float(div + 1) << "," << a * j / float(div + 1) - x << "," << l << "\n"
				<< 0.0f << "," << 0.0f << "," << 1.0f << "\n"
				<< 0.0f << "," << 0.0f << "\n"
				<< -c + (2 * c * (i + 1)) / float(div + 1) << "," << a * (j + 1) / float(div + 1) - x << "," << l << "\n"
				<< 0.0f << "," << 0.0f << "," << 1.0f << "\n"
				<< 1.0f << "," << 1.0f << "\n"
				<< -c + (2 * c * i) / float(div + 1) << "," << a * (j + 1) / float(div + 1) - x << "," << l << "\n"
				<< 0.0f << "," << 0.0f << "," << 1.0f << "\n"
				<< 0.0f << "," << 1.0f << "\n"

				//lado de trás 
				<< -c + (2 * c * i) / float(div + 1) << "," << a * j / float(div + 1) - x << "," << -l << "\n"
				<< 0.0f << "," << 0.0f << "," << -1.0f << "\n"
				<< 1.0f << "," << 0.0f << "\n"
				<< -c + (2 * c * (i + 1)) / float(div + 1) << "," << a * (j + 1) / float(div + 1) - x << "," << -l << "\n"
				<< 0.0f << "," << 0.0f << "," << -1.0f << "\n"
				<< 0.0f << "," << 1.0f << "\n"
				<< -c + (2 * c * (i + 1)) / float(div + 1) << "," << a * j / float(div + 1) - x << "," << -l << "\n"
				<< 0.0f << "," << 0.0f << "," << -1.0f << "\n"
				<< 0.0f << "," << 0.0f << "\n"

				<< -c + (2 * c * i) / float(div + 1) << "," << a * j / float(div + 1) - x << "," << -l << "\n"
				<< 0.0f << "," << 0.0f << "," << -1.0f << "\n"
				<< 1.0f << "," << 0.0f << "\n"
				<< -c + (2 * c * i) / float(div + 1) << "," << a * (j + 1) / float(div + 1) - x << "," << -l << "\n"
				<< 0.0f << "," << 0.0f << "," << -1.0f << "\n"
				<< 1.0f << "," << 1.0f << "\n"
				<< -c + (2 * c * (i + 1)) / float(div + 1) << "," << a * (j + 1) / float(div + 1) - x << "," << -l << "\n"
				<< 0.0f << "," << 0.0f << "," << -1.0f << "\n"
				<< 0.0f << "," << 1.0f << "\n"

				//lado da esquerda
				<< -c << "," << a * j / float(div + 1) - x << "," << -l + (2 * l * i) / float(div + 1) << "\n"
				<< -1.0f << "," << 0.0f << "," << 0.0f << "\n"
				<< 0.0f << "," << 1.0f << "\n"
				<< -c << "," << a * j / float(div + 1) - x << "," << -l + (2 * l * (i + 1)) / float(div + 1) << "\n"
				<< -1.0f << "," << 0.0f << "," << 0.0f << "\n"
				<< 0.0f << "," << 0.0f << "\n"
				<< -c << "," << a * (j + 1) / float(div + 1) - x << "," << -l + (2 * l * (i + 1)) / float(div + 1) << "\n"
				<< -1.0f << "," << 0.0f << "," << 0.0f << "\n"
				<< 1.0f << "," << 0.0f << "\n"

				<< -c << "," << a * j / float(div + 1) - x << "," << -l + (2 * l * i) / float(div + 1) << "\n"
				<< -1.0f << "," << 0.0f << "," << 0.0f << "\n"
				<< 0.0f << "," << 1.0f << "\n"
				<< -c << "," << a * (j + 1) / float(div + 1) - x << "," << -l + (2 * l * (i + 1)) / float(div + 1) << "\n"
				<< -1.0f << "," << 0.0f << "," << 0.0f << "\n"
				<< 1.0f << "," << 0.0f << "\n"
				<< -c << "," << a * (j + 1) / float(div + 1) - x << "," << -l + (2 * l * i) / float(div + 1) << "\n"
				<< -1.0f << "," << 0.0f << "," << 0.0f << "\n"
				<< 1.0f << "," << 1.0f << "\n"

				//lado da direita 
				<< c << "," << a * j / float(div + 1) - x << "," << -l + (2 * l * i) / float(div + 1) << "\n"
				<< 1.0f << "," << 0.0f << "," << 0.0f << "\n"
				<< 1.0f << "," << 1.0f << "\n"
				<< c << "," << a * (j + 1) / float(div + 1) - x << "," << -l + (2 * l * (i + 1)) / float(div + 1) << "\n"
				<< 1.0f << "," << 0.0f << "," << 0.0f << "\n"
				<< 0.0f << "," << 0.0f << "\n"
				<< c << "," << a * j / float(div + 1) - x << "," << -l + (2 * l * (i + 1)) / float(div + 1) << "\n"
				<< 1.0f << "," << 0.0f << "," << 0.0f << "\n"
				<< 1.0f << "," << 0.0f << "\n"

				<< c << "," << a * j / float(div + 1) - x << "," << -l + (2 * l * i) / float(div + 1) << "\n"
				<< 1.0f << "," << 0.0f << "," << 0.0f << "\n"
				<< 1.0f << "," << 1.0f << "\n"
				<< c << "," << a * (j + 1) / float(div + 1) - x << "," << -l + (2 * l * i) / float(div + 1) << "\n"
				<< 1.0f << "," << 0.0f << "," << 0.0f << "\n"
				<< 0.0f << "," << 1.0f << "\n"
				<< c << "," << a * (j + 1) / float(div + 1) - x << "," << -l + (2 * l * (i + 1)) / float(div + 1) << "\n"
				<< 1.0f << "," << 0.0f << "," << 0.0f << "\n"
				<< 0.0f << "," << 0.0f << "\n"
				;
		}
	}
	f.close();
}

vector<float> third_degree_bezier(float t, vector<vector<float>> cp){

	float p01x = (1 - t) * cp[0][0] + t * cp[1][0];
	float p01y = (1 - t) * cp[0][1] + t * cp[1][1];
	float p01z = (1 - t) * cp[0][2] + t * cp[1][2];

	float p12x = (1 - t) * cp[1][0] + t * cp[2][0];
	float p12y = (1 - t) * cp[1][1] + t * cp[2][1];
	float p12z = (1 - t) * cp[1][2] + t * cp[2][2];

	float p23x = (1 - t) * cp[2][0] + t * cp[3][0];
	float p23y = (1 - t) * cp[2][1] + t * cp[3][1];
	float p23z = (1 - t) * cp[2][2] + t * cp[3][2];

	float p012x = (1 - t) * p01x + t * p12x;
	float p012y = (1 - t) * p01y + t * p12y;
	float p012z = (1 - t) * p01z + t * p12z;

	float p123x = (1 - t) * p12x + t * p23x;
	float p123y = (1 - t) * p12y + t * p23y;
	float p123z = (1 - t) * p12z + t * p23z;

	vector<float> p;
	p.push_back(((1 - t) * p012x + t * p123x));
	p.push_back(((1 - t) * p012y + t * p123y));
	p.push_back(((1 - t) * p012z + t * p123z));

	return p;
}

vector<float> bezierPoint(float u, float v, vector<vector<float>> cp0, vector<vector<float>> cp1, vector<vector<float>> cp2, vector<vector<float>> cp3){

	vector<vector<float>> P;

	P.push_back(third_degree_bezier(u, cp0));
	P.push_back(third_degree_bezier(u, cp1));
	P.push_back(third_degree_bezier(u, cp2));
	P.push_back(third_degree_bezier(u, cp3));

	return third_degree_bezier(v, P);
}

void beziernormals(const char* file) {
	ifstream fin;
	string texto;
	fin.open(file);
	vector<float> points;
	while (getline(fin, texto)) {
		points = splitfloat(texto, ',', points);
	}
	for (int i = 0; i < points.size(); i++) {
		if (i == 0) continue;
	}
}

void bezier(float tesselacao, const char* bezier_file, const char* file){
	ifstream fin;
	string texto;
	vector<vector<float>> patches;
	vector<vector<float>> control_points;
	vector<float> indices;
	fin.open(bezier_file);

	getline(fin,texto);
	int number_patches = stoi(texto);

	for (int i = 0; i < number_patches; i++){
		getline(fin,texto);
		indices = splitfloat(texto, ',', indices);
		patches.push_back(indices);
		indices.clear();
	}

	getline(fin,texto);
	int number_control_points = stoi(texto);

	for (int i = 0; i < number_control_points; i++){
		getline(fin,texto);
		indices = splitfloat(texto, ',', indices);
		control_points.push_back(indices);
		indices.clear();
	}
	fin.close();

	vector<vector<float>> pontos;
	float u1, u2, v1, v2, incremento = 1.0f/tesselacao;
	for (int i = 0; i < number_patches; i++){
		for (int u = 0; u < tesselacao; u++){
			for (int v = 0; v < tesselacao; v++){
				u1 = u * incremento;
		        v1 = v * incremento;
		        u2 = (u + 1) * incremento;
		        v2 = (v + 1) * incremento;

				vector<vector<float>> cp0{  control_points[patches[i][0]],
										    control_points[patches[i][1]],
											control_points[patches[i][2]],
											control_points[patches[i][3]]
										};
				vector<vector<float>> cp1{ 
											control_points[patches[i][4]],
											control_points[patches[i][5]],
											control_points[patches[i][6]],
											control_points[patches[i][7]]
										};
				vector<vector<float>> cp2{  control_points[patches[i][8]],
											control_points[patches[i][9]],
											control_points[patches[i][10]],
											control_points[patches[i][11]]
										};
				vector<vector<float>> cp3{
											control_points[patches[i][12]],
											control_points[patches[i][13]],
											control_points[patches[i][14]],
											control_points[patches[i][15]]
										};

				vector<float> A = bezierPoint(u1, v1, cp0, cp1, cp2, cp3);
				vector<float> B = bezierPoint(u2, v1, cp0, cp1, cp2, cp3);
				vector<float> C = bezierPoint(u1, v2, cp0, cp1, cp2, cp3);
				vector<float> D = bezierPoint(u2, v2, cp0, cp1, cp2, cp3);

				vector<float> v1{ D[0] - A[0], D[1] - A[1], D[2] - A[2] };
				vector<float> v2{ C[0] - B[0], C[1] - B[1], C[2] - B[2] };
				vector<float> norm{ 0.0f, 0.0f, 0.0f };
				norm = cross(v1, v2, norm);
				norm = normalize(norm);

				pontos.push_back(A);
				pontos.push_back(norm);
				pontos.push_back({ 0.0f, 1.0f });
				pontos.push_back(B);
				pontos.push_back(norm);
				pontos.push_back({ 0.0f, 0.0f });
				pontos.push_back(D);
				pontos.push_back(norm);
				pontos.push_back({ 1.0f, 0.0f });

				pontos.push_back(A);
				pontos.push_back(norm);
				pontos.push_back({ 0.0f, 1.0f });
				pontos.push_back(D);
				pontos.push_back(norm);
				pontos.push_back({ 1.0f, 0.0f });	
				pontos.push_back(C);
				pontos.push_back(norm);
				pontos.push_back({ 1.0f, 1.0f });

				cp0.clear();
				cp1.clear();
				cp2.clear();
				cp3.clear();
			}
		}
	}
	ofstream fout;
	fout.open(file);
	for(auto i: pontos){
		int x = 0;
		while (x < i.size()) {
			fout << i[x];
			x++;
			if (x == i.size()) fout << "\n";
			else fout << ",";
		}
	}
	fout.close();
}

int main(int argc, char** argv) {
	if (argc <= 1) printf("Argumentos Insuficientes\n");
	string pasta = "teste/";
	string figura = argv[1];
	const char* caminho;
	if(figura != "bezier" && figura != "plano" && figura != "caixa" && figura != "esfera" && figura != "cilindro" && figura != "cone" && figura != "torus"){
		printf("Nome da figura invalido: %s\n", argv[1]);
	}
	else {
		if(figura == "bezier"){
			if(argc != 5){
				printf("Numero errado de argumentos\n");
			}
			else{
				pasta.append(argv[4]);
				string pasta1="teste/";
				pasta1.append(argv[3]);
				caminho = pasta.c_str();
				bezier(stof(argv[2]), pasta1.c_str(), caminho);
			}
		}
		if (figura == "torus"){
			if(argc != 7){
				printf("Numero errado de argumentos\n");
			}
			else{
				pasta.append(argv[6]);
				caminho = pasta.c_str();
				torus(stof(argv[2]), stof(argv[3]), stof(argv[4]), stof(argv[5]), caminho);
			}
		}
		if (figura == "plano") {
			if (argc != 5) printf("Numero errado de argumentos\n");
			else {
				pasta.append(argv[4]);
				caminho = pasta.c_str();
				plano(stof(argv[2]), stof(argv[3]), caminho);
			}
		}
		if (figura == "caixa") {
			if (argc != 7) printf("Numero errado de argumentos\n");
			else {
				pasta.append(argv[6]);
				caminho = pasta.c_str();
				caixa(stof(argv[2]), stof(argv[3]), stof(argv[4]), stof(argv[5]), caminho);
			}
		}
		if (figura == "esfera") {
			if (argc != 6) printf("Numero errado de argumentos\n");
			else {
				pasta.append(argv[5]);
				caminho = pasta.c_str();
				esfera(stof(argv[2]), stof(argv[3]), stof(argv[4]), caminho);
			}
		}
		if (figura == "cone") {
			if (argc != 7) printf("Numero errado de argumentos\n");
			else {
				pasta.append(argv[6]);
				caminho = pasta.c_str();
				cone(stof(argv[2]), stof(argv[3]), stof(argv[4]), stof(argv[5]), caminho);
			}
		}
		if (figura == "cilindro") {
			if (argc != 6) printf("Numero errado de argumentos\n");
			else {
				pasta.append(argv[5]);
				caminho = pasta.c_str();
				cilindro(stof(argv[2]), stof(argv[3]), stoi(argv[4]), caminho);
			}
		}
	}
}