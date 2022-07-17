#include <windows.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#include <IL/il.h>
#endif
#define _USE_MATH_DEFINES
#include "headers/tinyxml2.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>


using namespace std;
using namespace tinyxml2;

unsigned int picked = 0;
bool pressing_left_button = false, pressing_right_button = false, break_xml = false, normals = false, fpscamera = false;
vector<string> files;
vector<vector<float>> lightpos, drawnormals;
int timebase, ind = 0, tracking = 0, w, h;
int caixa = -1, cilindro = -1, cone = -1, esfera = -1, plano = -1, torus = -1, bezier = -1,
toruspicking = -1, teapot = -1, sol = -1, mercurio = -1, venus = -1, terra = -1, lua = -1, marte = -1, jupiter = -1, saturno = -1, urano = -1, neptuno = -1, plutao = -1;
GLenum lightn[8] = { GL_LIGHT0, GL_LIGHT1, GL_LIGHT2, GL_LIGHT3, GL_LIGHT4, GL_LIGHT5, GL_LIGHT6, GL_LIGHT7 };
float frames, camX = 00, camY = 30, camZ = 40, Px = 0, Py = 0, Pz = 0;
GLuint vertexbuffers[10], normalbuffers[10], texturebuffers[10];
int verticeCount[512], triangleCount = 0;
int startX, startY, alpha = 0, beta = 45, r = 50;
char text[64];
string xml;
float global_time = 1;

//Retorna o planeta que foi selecionado
string qual_objeto_picking(int pick) {
	if (pick == sol) return "Sol";
	else if (pick == mercurio) return "Mercurio";
	else if (pick == venus) return "Venus";
	else if (pick == terra) return "Terra";
	else if (pick == lua) return "Lua";
	else if (pick == marte) return "Marte";
	else if (pick == jupiter) return "Jupiter";
	else if (pick == saturno) return "Saturno";
	else if (pick == urano) return "Urano";
	else if (pick == neptuno) return "Neptuno";
	else if (pick == plutao) return "Plutao";
	else if (pick == toruspicking) return "Torus";
	else if (pick == teapot) return "Bezier";
	else return "";
}

//Retorna o buffer associado ao objeto
int qual_objeto(string object) {
	if (object == "teste/caixa.3d") return caixa;
	if (object == "teste/cilindro.3d") return cilindro;
	if (object == "teste/cone.3d") return cone;
	if (object == "teste/esfera.3d") return esfera;
	if (object == "teste/plano.3d") return plano;
	if (object == "teste/torus.3d") return torus;
	if (object == "teste/bezier.3d") return bezier;
	return -1;
}

void buildRotMatrix(float* x, float* y, float* z, float* m) {
	m[0] = x[0]; m[1] = x[1]; m[2] = x[2]; m[3] = 0;
	m[4] = y[0]; m[5] = y[1]; m[6] = y[2]; m[7] = 0;
	m[8] = z[0]; m[9] = z[1]; m[10] = z[2]; m[11] = 0;
	m[12] = 0; m[13] = 0; m[14] = 0; m[15] = 1;
}

void cross(float* a, float* b, float* res) {

	res[0] = a[1] * b[2] - a[2] * b[1];
	res[1] = a[2] * b[0] - a[0] * b[2];
	res[2] = a[0] * b[1] - a[1] * b[0];
}

void normalize(float* a) {

	float l = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
	a[0] = a[0] / l;
	a[1] = a[1] / l;
	a[2] = a[2] / l;
}

void multMatrixVector(float* m, float* v, float* res) {

	for (int j = 0; j < 4; ++j) {
		res[j] = 0;
		for (int k = 0; k < 4; ++k) {
			res[j] += v[k] * m[j * 4 + k];
		}
	}
}

void getCatmullRomPoint(float t, vector<float> p0, vector<float> p1, vector<float> p2, vector<float> p3, float* pos, float* deriv) {
	// catmull-rom matrix
	float M[4][4] = { {-0.5f,  1.5f, -1.5f,  0.5f},
						{ 1.0f, -2.5f,  2.0f, -0.5f},
						{-0.5f,  0.0f,  0.5f,  0.0f},
						{ 0.0f,  1.0f,  0.0f,  0.0f}
	};

	float T[4] = { pow(t,3.0f), pow(t,2.0f), t, 1.0f };
	float Td[4] = { 3.0f * pow(t,2.0f), 2.0f * t, 1.0f, 0.0f };

	float X[4] = { p0[0], p1[0], p2[0], p3[0] };
	float Y[4] = { p0[1], p1[1], p2[1], p3[1] };
	float Z[4] = { p0[2], p1[2], p2[2], p3[2] };

	// Compute A = M * P
	float A[4][4] = { 0 };

	multMatrixVector(*M, X, A[0]);
	multMatrixVector(*M, Y, A[1]);
	multMatrixVector(*M, Z, A[2]);

	// Compute pos = A * T
	multMatrixVector(*A, T, pos);

	// compute deriv = A * Td
	multMatrixVector(*A, Td, deriv);
}

// given  global t, returns the point in the curve
void getGlobalCatmullRomPoint(float gt, float* pos, float* deriv, vector<vector<float>> pontos) {

	int point_count = pontos.size();
	float t = gt * point_count; // this is the real global t
	int index = floor(t);  // which segment
	t = t - index; // where within  the segment

	// indices store the points
	int indices[4];
	indices[0] = (index + point_count - 1) % point_count;
	indices[1] = (indices[0] + 1) % point_count;
	indices[2] = (indices[1] + 1) % point_count;
	indices[3] = (indices[2] + 1) % point_count;
	getCatmullRomPoint(t, pontos[indices[0]], pontos[indices[1]], pontos[indices[2]], pontos[indices[3]], pos, deriv);
}

//Desenha as curvas Catmull-Rom
void renderCatmullRomCurve(float tesselacao, vector<vector<float>> pontos) {
	float pos[4];
	float deriv[4];
	glDisable(GL_LIGHTING);
	glBegin(GL_LINE_LOOP);
	glColor3f(0.7f, 0.7f, 0.7f);
	for (int t = 0; t < tesselacao; t++) {
		getGlobalCatmullRomPoint(float(t / 60.0), pos, deriv, pontos);
		glVertex3f(pos[0], pos[1], pos[2]);
	}
	glEnd();
	glEnable(GL_LIGHTING);
}

//Retorna um vetor de floats com os valores da string s separados pelo delimitador d
vector<float> splitfloat(const string& s, char d, vector<float> vetor) {
	stringstream ss(s);
	string texto;
	while (getline(ss, texto, d)) {
		vetor.push_back(stof(texto));
	}
	return vetor;
}

//Retorna um vetor de strings com os valores da string s separados pelo delimitador d
vector<string> splitstring(const string& s, char d, vector<string> vetor) {
	stringstream ss(s);
	string texto;
	while (getline(ss, texto, d)) {
		vetor.push_back(texto);
	}
	return vetor;
}

//Faz a translação 
void translacao(string instruction) {
	float pos[4], deriv[4], z[4], y[4] = { 0,1,0 }, m[16];
	vector<float> x;
	x = splitfloat(&instruction[2], ',', x);
	float time = x[0];
	x.erase(x.begin());

	vector<vector<float>> pontos;
	vector<float> pon;
	for (int i = 0; i < x.size(); i = i + 3) {
		pon.push_back(x[i]);
		pon.push_back(x[i + 1]);
		pon.push_back(x[i + 2]);
		pontos.push_back(pon);
		pon.clear();
	}

	renderCatmullRomCurve(60, pontos);

	float current_time = glutGet(GLUT_ELAPSED_TIME) / 1000.0;
	float tempo = current_time * global_time / time;
	getGlobalCatmullRomPoint(tempo, pos, deriv, pontos);

	normalize(deriv);
	cross(deriv, y, z);
	normalize(z);
	cross(z, deriv, y);
	normalize(y);

	buildRotMatrix(deriv, y, z, m);
	glPushMatrix();
	glTranslatef(pos[0], pos[1], pos[2]);
	glMultMatrixf(m);
}

//Faz a rotação
void rotacao(vector<float> arg) {
	float time = arg[0];
	float current_time = glutGet(GLUT_ELAPSED_TIME) / 1000.0;
	float tempo = 360.0f * current_time / time;
	tempo = tempo * global_time;
	glPushMatrix();
	glRotatef(tempo, arg[1], arg[2], arg[3]);
}

//Desenha os planetas com uma cor sólida cujo valor lhes é atribuído
void desenhaPicking(int code, string model) {
	glDisable(GL_LIGHTING);
	float color = code / 255.0f;
	int object = qual_objeto(model);
	glColor3f(color, color, color);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffers[object]);
	glVertexPointer(3, GL_FLOAT, 0, 0);
	glDrawArrays(GL_TRIANGLES, 0, verticeCount[object]);

	if (model == "teste/esfera.3d") { 
		if (sol == -1) sol = code;
		else if (mercurio == -1) mercurio = code;
		else if (venus == -1) venus = code;
		else if (terra == -1) terra = code;
		else if (lua == -1) lua = code;
		else if (marte == -1) marte = code;
		else if (jupiter == -1) jupiter = code;
		else if (saturno == -1) saturno = code;
		else if (urano == -1) urano = code;
		else if (neptuno == -1) neptuno = code;
		else if (plutao == -1) plutao = code;
	}
	else if (model == "teste/bezier.3d") teapot = code;
	else if (model == "teste/torus.3d") toruspicking = code;
}

//Executa uma instrução
void executa_instrucao(string instruction) {
	vector<float> arg;
	string texto;
	vector<string> argf;
	switch (instruction[0]) {
	case('T'):
		if (instruction == "T,inv") {
			glPopMatrix();
		}
		else {
			translacao(instruction);
		}
		break;
	case('R'):
		if (instruction[2] == 'A') {
			arg = splitfloat(&instruction[4], ',', arg);
			glRotatef(arg[0], arg[1], arg[2], arg[3]);
		}
		else {
			if (instruction[2] == 'T') {
				arg = splitfloat(&instruction[4], ',', arg);
				rotacao(arg);
			}
			else {
				glPopMatrix();
			}
		}
		break;
	case('S'):
		arg = splitfloat(&instruction[2], ',', arg);
		glScalef(arg[0], arg[1], arg[2]);
		break;
	case('C'):
		if (instruction[1] == ',') {
			arg = splitfloat(&instruction[2], ',', arg);
			glColor3f(arg[0], arg[1], arg[2]);
		}
		else {
			glColor3f(1.0f, 1.0f, 1.0f);
		}
		break;
	default:
		if(instruction.find("teste") == 0){
			argf = splitstring(instruction, ',', argf);
			int object = qual_objeto(argf[0]);
			float inc = 0;
			if (argf[1] == "T") inc = 1;

			float diff[4] = { stof(argf[2 + inc]), stof(argf[3 + inc]), stof(argf[4 + inc]), 1.0 };
			float spec[4] = { stof(argf[5 + inc]), stof(argf[6 + inc]), stof(argf[7 + inc]), 1.0 };
			float emi[4] = { stof(argf[8 + inc]), stof(argf[9 + inc]), stof(argf[10 + inc]), 1.0 };
			float amb[4] = { stof(argf[11 + inc]), stof(argf[12 + inc]), stof(argf[13 + inc]), 1.0 };

			glMaterialfv(GL_FRONT, GL_DIFFUSE, diff);
			glMaterialfv(GL_FRONT, GL_SPECULAR, spec);
			glMaterialfv(GL_FRONT, GL_EMISSION, emi);
			glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
			glMaterialf(GL_FRONT, GL_SHININESS, 0);

			GLuint id;
			if (argf[1] == "T") id = stoi(argf[2]);
			else id = 0;

			glBindBuffer(GL_ARRAY_BUFFER, vertexbuffers[object]);
			glVertexPointer(3, GL_FLOAT, 0, 0);

			glBindBuffer(GL_ARRAY_BUFFER, normalbuffers[object]);
			glNormalPointer(GL_FLOAT, 0, 0);

			glBindBuffer(GL_ARRAY_BUFFER, texturebuffers[object]);
			glTexCoordPointer(2, GL_FLOAT, 0, 0);

			glBindTexture(GL_TEXTURE_2D, id);
			glDrawArrays(GL_TRIANGLES, 0, verticeCount[object]/3);
			glBindTexture(GL_TEXTURE_2D, 0);

			if (normals) {
				glDisable(GL_LIGHTING);
				for (int i = 0; i < drawnormals[object].size(); i += 6) {
					glColor3f(0.0f, 0.8f, 0.0f); 
					glBegin(GL_LINES);
					glVertex3f(drawnormals[object][i], drawnormals[object][i + 1], drawnormals[object][i + 2]);
					glVertex3f(drawnormals[object][i + 3], drawnormals[object][i + 4], drawnormals[object][i + 5]);
					glEnd();
				}
				glEnable(GL_LIGHTING);
			}
		}
		break;
	}
}

//Calcula o valor da componente vermelha do planeta slecionado
unsigned char picking(int x, int y) {
	vector<string> argf;
	GLint viewport[4];
	unsigned char res[4];
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity();
	if (fpscamera) {
		gluLookAt(Pz, Py, Px,
			camZ, camY, camX,
			0.0f, 1.0f, 0.0f);
	}
	else {
		gluLookAt(camZ, camY, camX,
			0.0, 0.0, 0.0,
			0.0f, 1.0f, 0.0f);
	}
	glDepthFunc(GL_LEQUAL);
	for (int i = 0; i < files.size(); i++) {
		if (files[i].find("teste") == 0) {
			argf = splitstring(files[i], ',', argf);
			desenhaPicking(i + 1, argf[0]);
			argf.clear();
		}
		else {
			executa_instrucao(files[i]);
		}
	}
	glDepthFunc(GL_LESS);
	glGetIntegerv(GL_VIEWPORT, viewport);
	glReadPixels(x, viewport[3] - y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, res);
	glEnable(GL_LIGHTING);	
	glEnable(GL_TEXTURE_2D);
	return res[0];
}

void renderTexto() {
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0, w, h, 0);
	glMatrixMode(GL_MODELVIEW);
	glDisable(GL_DEPTH_TEST);
	glPushMatrix();
	glLoadIdentity();
	glRasterPos2d(10, 20);
	for (int i = 0; text[i] != '\0'; i++) {
		glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, text[i]);
	}
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	glEnable(GL_DEPTH_TEST);
}

void carregar_rato(int button, int state, int xx, int yy) {
	if (state == GLUT_DOWN) {
		startX = xx;
		startY = yy;
		string c;
		int i;
		if (button == GLUT_LEFT_BUTTON) tracking = 1;
		else if (button == GLUT_RIGHT_BUTTON) tracking = 2;
		else if (button == GLUT_MIDDLE_BUTTON && xml == "sistemasolar.xml") {
				toruspicking = -1;
				teapot = -1;
				sol = -1;
				mercurio = -1;
				venus = -1;
				terra = -1;
				lua = -1;
				marte = -1;
				jupiter = -1;
				saturno = -1;
				urano = -1;
				neptuno = -1;
				plutao = -1;
				picked = picking(xx, yy);
				if (picked) {
					c = "Escolhido: " + qual_objeto_picking(picked);
					i = 0;
					while(c[i] != '\0') {
						text[i] = c[i];
						i++;
					}
					text[i] = '\0';
				}
				else {
					c = "Nada selecionado";
					i = 0;
					while (c[i] != '\0') {
						text[i] = c[i];
						i++;
					}
					text[i] = '\0';
				}
				glutPostRedisplay();

			}
		else tracking = 0;
	}
	else if (state == GLUT_UP) {
		if (tracking == 1) {
			alpha += (xx - startX);
			if (fpscamera) beta += (startY - yy);
			else beta += (yy - startY);
		}
		else if (tracking == 2) {

			r -= yy - startY;
			if (r < 3)
				r = 3.0;
		}
		tracking = 0;
	}
}

void rato(int xx, int yy) {
	int deltaX, deltaY;
	int alphaAux, betaAux;
	int rAux;
	int y;

	if (!tracking) return;
	deltaX = xx - startX;
	if (fpscamera) deltaY = startY - yy;
	else deltaY = yy - startY;

	if (tracking == 1) {
		alphaAux = alpha + deltaX;
		betaAux = beta + deltaY;
		if (betaAux > 85.0)
			betaAux = 85.0;
		else if (betaAux < -85.0)
			betaAux = -85.0;

		rAux = r;
	}
	else if (tracking == 2) {
			alphaAux = alpha;
			betaAux = beta;
			if (!fpscamera) {
				rAux = r - deltaY;
				if (rAux < 3) rAux = 3;
			}
		}
	if (fpscamera) {
		camX = Px + sin(alphaAux * 3.14 / 180.0);
		camY = Py + sin(betaAux * 3.14 / 180.0);
		camZ = Pz + cos(alphaAux * 3.14 / 180.0);
	}
	else {
		camX = rAux * sin(alphaAux * 3.14 / 180.0) * cos(betaAux * 3.14 / 180.0);
		camZ = rAux * cos(alphaAux * 3.14 / 180.0) * cos(betaAux * 3.14 / 180.0);
		camY = rAux * sin(betaAux * 3.14 / 180.0);
	}
}

void input(unsigned char key, int x, int y) {
	float dx = camX - Px, dz = camZ - Pz;
	switch (key) {
	case('o'):
		r += 0.1f;
		break;
	case('i'):
		if (r >= 0) r -= 0.1f;
		break;
	case 'l':
		glPolygonMode(GL_FRONT, GL_LINE);
		break;
	case 'p':
		glPolygonMode(GL_FRONT, GL_POINT);
		break;
	case 'f':
		glPolygonMode(GL_FRONT, GL_FILL);
		break;
	case 'n':
		if (normals) normals = false;
		else normals = true;
		break;
	case('c'):
		if (fpscamera) {
			fpscamera = false;
			camX = Px;
			camY = Py;
			camZ = Pz;
		}
		else {
			fpscamera = true;
			Px = camX;
			Py = camY;
			Pz = camZ;
			camX = 0;
			camY = 0;
			camZ = 0;
		}
		break;
	case('e'):
		if (fpscamera) {
			Py += 1;
			camY += 1;
		}
		break;
	case('q'):
		if (fpscamera) {
			Py += -1;
			camY += -1;
		}
		break;
	case('w'):          
		if (fpscamera) {
			Px += dx * 1.5;
			Pz += dz * 1.5;
			camX += dx * 1.5;
			camZ += dz * 1.5;
		}
		break;
	case('s'):
		if (fpscamera) {
			Px -= dx * 1.5;
			Pz -= dz * 1.5;
			camX -= dx * 1.5;
			camZ -= dz * 1.5;
		}
		break;
	case('a'):
		if (fpscamera) {
			Px -= dz * 1.5;
			Pz += dx * 1.5;
			camX -= dz * 1.5;
			camZ -= -dx * 1.5;
		}
		break;
	case('d'):
		if (fpscamera) {
			Px += dz * 1.5;
			Pz -= dx * 1.5;
			camX += dz * 1.5;
			camZ += -dx * 1.5;
		}
		break;
	case('x'):
		if (global_time < 20.0f) global_time += 0.2f;
		break;
	case('z'):
		if (global_time > 0.2f) global_time -= 0.2f;
		break;
	default:
		break;
	}
	glutPostRedisplay();
}

//Preenche os buffers
void escreve(const char* file) {
	verticeCount[ind] = 0;
	ifstream f;
	string texto;
	string vnt = "vertex";
	vector<float> v, n, dn, t;
	f.open(file);

	while (getline(f, texto)) {
		if (vnt == "vertex") {
			v = splitfloat(texto, ',', v);
			vnt = "normal";
		}
		else {
			if (vnt == "normal") {
				n = splitfloat(texto, ',', n);
				vnt = "texture";
				
			}
			else {
				t = splitfloat(texto, ',', t);
				vnt = "vertex";
			}
		}
	}
	f.close();
	verticeCount[ind] = v.size();
	for (int i = 0; i < v.size(); i += 3) {
		float ponto[3] = { v[i], v[i + 1], v[i + 2] };
		float normal[3] = { n[i], n[i + 1], n[i + 2] };
		dn.push_back(v[i]);
		dn.push_back(v[i + 1]);
		dn.push_back(v[i + 2]);
		dn.push_back(v[i] + n[i]);
		dn.push_back(v[i + 1] + n[i + 1]);
		dn.push_back(v[i + 2] + n[i + 2]);
	}
	drawnormals.push_back(dn);
	dn.clear();
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffers[ind]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * v.size(), v.data(), GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, normalbuffers[ind]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * n.size(), n.data(), GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, texturebuffers[ind]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * t.size(), t.data(), GL_STATIC_DRAW);
}

//Abre um caminho especificado para leitura
void desenhaaux(string path) {
	ifstream f;
	f.open(path);
	string diretorias;
	getline(f, diretorias);
	escreve(path.c_str());
	f.close();
}

//Verifica se o modelo em questão já foi lido antes. Se não foi lido, atribui um inteiro ao modelo correspondente a um buffer
void desenha(vector<string> files) {
	vector<string> argf;
	for (int i = 0; i < files.size(); i++) {
		if (files[i][1] != ',' && files[i] != "Cor_predefinida") {
			argf = splitstring(files[i], ',', argf);
			if (caixa == -1 && argf[0] == "teste/caixa.3d") {
				desenhaaux("teste/caixa.3d");
				caixa = ind;
				ind++;
			}
			if (cilindro == -1 && argf[0] == "teste/cilindro.3d") {
				desenhaaux("teste/cilindro.3d");
				cilindro = ind;
				ind++;
			}
			if (cone == -1 && argf[0] == "teste/cone.3d") {
				desenhaaux("teste/cone.3d");
				cone = ind;
				ind++;
			}
			if (esfera == -1 && argf[0] == "teste/esfera.3d") {
				desenhaaux("teste/esfera.3d");
				esfera = ind;
				ind++;
			}
			if (plano == -1 && argf[0] == "teste/plano.3d") {
				desenhaaux("teste/plano.3d");
				plano = ind;
				ind++;
			}
			if (torus == -1 && argf[0] == "teste/torus.3d") {
				desenhaaux("teste/torus.3d");
				torus = ind;
				ind++;
			}
			if (bezier == -1 && argf[0] == "teste/bezier.3d") {
				desenhaaux("teste/bezier.3d");
				bezier = ind;
				ind++;
			}
			argf.clear();
		}
	}
}

void changeSize(int ww, int hh) {
	// Prevent a divide by zero, when window is too short
	// (you can’t make a window with zero width).
	w = ww;
	h = hh;
	if (hh == 0)
		hh = 1;
	// compute window's aspect ratio
	float ratio = ww * 1.0f / hh;
	// Set the projection matrix as current
	glMatrixMode(GL_PROJECTION);
	// Load the identity matrix
	glLoadIdentity();
	// Set the viewport to be the entire window
	glViewport(0, 0, ww, hh);
	// Set the perspective
	gluPerspective(45.0f, ratio, 1.0f, 1000.0f);
	// return to the model view matrix mode
	glMatrixMode(GL_MODELVIEW);
}

//Mostra informações no título da janela e calcula os frames por segundo
void windowtitle() {
	frames++;
	int time = glutGet(GLUT_ELAPSED_TIME);
	if (time - timebase > 1000) {
		int fps = frames * 1000.0f / (time - timebase);
		char fpsstring[(((sizeof fps) * CHAR_BIT) + 2) / 3 + 50];
		sprintf(fpsstring, "%s    FPS: %d  Speed: x%.2f  Triangles: %d", xml.c_str(), fps, global_time, triangleCount);
		glutSetWindowTitle(fpsstring);
		timebase = time;
		frames = 0;
	}
}

//Inicializa e modifica as luzes a serem utilizadas
void lighinit() {
	vector<float> arg;
	vector<int> indtoremove;
	int l = 0;
	for (int i = 0; i < files.size(); i++) {
		if (files[i][0] == 'L') {
			indtoremove.push_back(i);
			glEnable(lightn[l]);
			arg = splitfloat(&files[i][4], ',', arg);
			GLfloat dark[4] = { 0.2, 0.2, 0.2, 1.0 };
			GLfloat white[4] = { 1.0, 1.0, 1.0, 1.0 };

			// light colors
			glLightfv(lightn[l], GL_AMBIENT, dark);
			glLightfv(lightn[l], GL_DIFFUSE, white);
			glLightfv(lightn[l], GL_SPECULAR, white);

			if (files[i][2] == 'P') {
				lightpos.push_back({ arg[0], arg[1], arg[2], 1.0 });
			}
			if (files[i][2] == 'D') {
				lightpos.push_back({ arg[0], arg[1], arg[2], 0.0 });
			}
			if (files[i][2] == 'S') {
				float dir[3] = { arg[4], arg[5], arg[6] };
				lightpos.push_back({arg[0], arg[1], arg[2], arg[3], arg[4], arg[5], arg[6] });
			}
			arg.clear();
		}
		l++;
	}
	int j = 0;
	for (auto i : indtoremove) {
		files.erase(files.begin() + i - j);
		j++;
	}
}

//Desenha eixos
void desenha_eixos() {
	glBegin(GL_LINES);
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(1000.0f, 0.0f, 0.0f);

	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 1000.0f, 0.0f);

	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 1000.0f);
	glEnd();
}

void renderScene(void) {

	// clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	float dark[] = { 0.2, 0.2, 0.2, 1.0 };
	float white[] = { 0.8, 0.8, 0.8, 1.0 };
	float red[] = { 0.8, 0.2, 0.2, 1.0 };

	// set camera
	glLoadIdentity();
	if (fpscamera) {
		gluLookAt(Pz, Py, Px,
			camZ, camY, camX,
			0.0f, 1.0f, 0.0f);
	}
	else {
		gluLookAt(camZ, camY, camX,
			0.0, 0.0, 0.0,
			0.0f, 1.0f, 0.0f);
	}
	for (int i = 0; i < lightpos.size(); i++) {
		if (lightpos[i].size() > 4) {
			float dir[3] = { lightpos[i][4], lightpos[i][5], lightpos[i][6] };
			float pos[4] = { lightpos[i][1], lightpos[i][2], lightpos[i][3], 1.0f };
			glLightf(lightn[i], GL_SPOT_CUTOFF, lightpos[i][0]);
			glLightfv(lightn[i], GL_SPOT_DIRECTION, dir);
			glLightfv(lightn[i], GL_POSITION, pos);
		}
		else {
			float pos[4] = { lightpos[i][0], lightpos[i][1], lightpos[i][2], lightpos[i][3] };
			glLightfv(lightn[i], GL_POSITION, pos);
		}
	}
	//desenha_eixos();
	for (int i = 0; i < files.size(); i++) {
		executa_instrucao(files[i]);
	}
	glColor3f(1.0f, 1.0f, 1.0f);
	glDisable(GL_LIGHTING);
	renderTexto();
	glEnable(GL_LIGHTING);
	windowtitle();

	// End of frame
	glutSwapBuffers();

}

//Faz parse do documento XML
void ler_group(XMLElement* elemento) {
	vector<string> inv;
	bool t = false, r = false, s = false, c = false, m = false;
	string group = "group", translate = "translate", rotate = "rotate", scale = "scale", models = "models", color = "color", point = "point";
	XMLElement* action = elemento->FirstChildElement();
	while (action) {
		string instrucao, instrucao_inv;
		if (action->Name() != point && action->Name() != translate && action->Name() != rotate && action->Name() != scale && action->Name() != color && action->Name() != models && action->Name() != group) {
			cout << "ERROR: UNKNOWN NODE \"" << action->Name() << "\"\n";
			action = action->NextSiblingElement();
			break_xml = true;
		}
		else {
			if (action->Name() == translate) {
				if (t) {
					break_xml = true;
					printf("ERROR: MULTIPLE TRANSLATE NODES IN A GROUP\n");
				}
				if (m) {
					break_xml = true;
					printf("ERROR: TRANSLATE NODE SHOULD BE IN THE BEGINNING OF THE GROUP\n");
				}
				XMLElement* point = action->FirstChildElement();
				instrucao = "T," + to_string(action->FloatAttribute("time"));
				int point_min = 0;
				while (point) {
					instrucao = instrucao + "," + to_string(point->FloatAttribute("X")) + ',' + to_string(point->FloatAttribute("Y")) + ',' + to_string(point->FloatAttribute("Z"));
					point = point->NextSiblingElement();
					point_min++;
				}
				if (point_min < 4) {
					printf("ERROR: TRANSLATE NOTE SHOULD HAVE 4 OR MORE POINTS.\n");
					break_xml = true;
				}
				instrucao_inv = "T,inv";
				files.push_back(instrucao);
				inv.push_back(instrucao_inv);
				t = true;
			}
			if (action->Name() == rotate) {
				if (r) {
					break_xml = true;
					printf("ERROR: MULTIPLE ROTATE NODES IN A GROUP\n");
				}
				if (m) {
					break_xml = true;
					printf("ERROR: ROTATE NODE SHOULD BE IN THE BEGINNING OF THE GROUP\n");
				}
				if (action->FloatAttribute("angle")) {
					instrucao = "R,A," + to_string(action->FloatAttribute("angle")) + ',' + to_string(action->FloatAttribute("axisX")) + ',' + to_string(action->FloatAttribute("axisY")) + ',' + to_string(action->FloatAttribute("axisZ"));
					instrucao_inv = "R,A," + to_string(-action->FloatAttribute("angle")) + ',' + to_string(action->FloatAttribute("axisX")) + ',' + to_string(action->FloatAttribute("axisY")) + ',' + to_string(action->FloatAttribute("axisZ"));
					files.push_back(instrucao);
					inv.push_back(instrucao_inv);
				}
				else {
					instrucao = "R,T," + to_string(action->FloatAttribute("time")) + ',' + to_string(action->FloatAttribute("axisX")) + ',' + to_string(action->FloatAttribute("axisY")) + ',' + to_string(action->FloatAttribute("axisZ"));
					instrucao_inv = "R,inv";
					files.push_back(instrucao);
					inv.push_back(instrucao_inv);
				}
				r = true;
			}
			if (action->Name() == scale) {
				if (s) {
					break_xml = true;
					printf("ERROR: MULTIPLE SCALE NODES IN A GROUP\n");
				}
				if (m) {
					break_xml = true;
					printf("ERROR: SCALE NODE SHOULD BE IN THE BEGINNING OF THE GROUP\n");
				}
				if (action->FloatAttribute("X")) instrucao = "S," + to_string(action->FloatAttribute("X")) + ",";
				else instrucao = "S,1,";
				if (action->FloatAttribute("Y")) instrucao = instrucao + to_string(action->FloatAttribute("Y")) + ",";
				else instrucao = instrucao + "1,";
				if (action->FloatAttribute("Z")) instrucao = instrucao + to_string(action->FloatAttribute("Z"));
				else instrucao = instrucao + "1";
				instrucao_inv = "S," + to_string(1 / action->FloatAttribute("X")) + ',' + to_string(1 / action->FloatAttribute("Y")) + ',' + to_string(1 / action->FloatAttribute("Z"));
				files.push_back(instrucao);
				inv.push_back(instrucao_inv);
				s = true;
			}
			if (action->Name() == color) {
				if (c) {
					break_xml = true;
					printf("ERROR: MULTIPLE COLOR NODES IN A GROUP\n");
				}
				if (m) {
					break_xml = true;
					printf("ERROR: COLOR NODE SHOULD BE IN THE BEGINNING OF THE GROUP\n");
				}
				instrucao = "C," + to_string(action->FloatAttribute("R")) + ',' + to_string(action->FloatAttribute("G")) + ',' + to_string(action->FloatAttribute("B"));
				instrucao_inv = "Cor_predefinida";
				files.push_back(instrucao);
				inv.push_back(instrucao_inv);
				c = true;
			}
			if (action->Name() == models) {
				if (m) {
					break_xml = true;
					printf("ERROR: MULTIPLE MODELS NODES IN A GROUP\n");
				}
				XMLElement* model = action->FirstChildElement();
				ifstream f;
				while (model) {
					instrucao = model->Attribute("file");
					f.open(instrucao);
					if (!f.is_open()) {
						cout << "Path \"" << instrucao << "\" was not found\n";
						break_xml = true;
					}
					else {
						if (model->Attribute("texture")) instrucao = instrucao + ",T," + model->Attribute("texture") + ",";
						else instrucao = instrucao + ",C,";
						if (model->Attribute("diffR")) instrucao = instrucao + to_string(model->FloatAttribute("diffR")) + ",";
						else instrucao = instrucao + "0.8,";
						if (model->Attribute("diffG")) instrucao = instrucao + to_string(model->FloatAttribute("diffG")) + ",";
						else instrucao = instrucao + "0.8,";
						if (model->Attribute("diffB")) instrucao = instrucao + to_string(model->FloatAttribute("diffB")) + ",";
						else instrucao = instrucao + "0.8,";
						
						instrucao = instrucao + to_string(model->FloatAttribute("specR")) + ',' + to_string(model->FloatAttribute("specG")) + ',' + to_string(model->FloatAttribute("specB")) + ',';
						instrucao = instrucao + to_string(model->FloatAttribute("emiR")) + ',' + to_string(model->FloatAttribute("emiG")) + ',' + to_string(model->FloatAttribute("emiB")) + ',';


						if (model->Attribute("ambR")) instrucao = instrucao + to_string(model->FloatAttribute("ambR")) + ",";
						else instrucao = instrucao + "0.2,";
						if (model->Attribute("ambG")) instrucao = instrucao + to_string(model->FloatAttribute("ambG")) + ",";
						else instrucao = instrucao + "0.2,";
						if (model->Attribute("ambB")) instrucao = instrucao + to_string(model->FloatAttribute("ambB"));
						else instrucao = instrucao + "0.2";

						files.push_back(instrucao); 
					}
					model = model->NextSiblingElement();
				}
				m = true;
			}
			if (action->Name() == group) {
				ler_group(action);
			}
			action = action->NextSiblingElement();
		}
	}
	for (int i = inv.size() - 1; i >= 0; i--) {
		files.push_back(inv[i]);
	}
}

//Faz parse da componente da luz do documento XML
void ler_lights(XMLElement* elemento) {
	string type, instrucao;
	XMLElement* light = elemento->FirstChildElement();
	while (light) {
		type = light->Attribute("type");
		if (type == "POINT") {
			instrucao = "L,P," + to_string(light->FloatAttribute("posX")) + "," + to_string(light->FloatAttribute("posY")) + "," + to_string(light->FloatAttribute("posZ"));
		}
		if (type == "DIRECTIONAL") {
			instrucao = "L,D," + to_string(light->FloatAttribute("vecX")) + "," + to_string(light->FloatAttribute("vecY")) + "," + to_string(light->FloatAttribute("vecZ"));
		}
		if (type == "SPOT") {
			instrucao = "L,S," + to_string(light->FloatAttribute("angle")) + ',' + to_string(light->FloatAttribute("posX")) + "," + to_string(light->FloatAttribute("posY")) + "," + to_string(light->FloatAttribute("posZ")) + ',' + to_string(light->FloatAttribute("vecX")) + "," + to_string(light->FloatAttribute("vecY")) + "," + to_string(light->FloatAttribute("vecZ"));
		}
		files.push_back(instrucao);
		light = light->NextSiblingElement();
	}
}

//Testa se o documento XML tem erros
void leitor_xml(string arg) {
	string pasta = "modelos/";
	string group = "group", lights = "lights";
	pasta.append(arg);
	char* file = new char[pasta.size() + 1];
	copy(pasta.begin(), pasta.end(), file);
	file[pasta.size()] = '\0';
	tinyxml2::XMLDocument doc;
	tinyxml2::XMLError flag = doc.LoadFile(file);
	if (!flag) {
		XMLElement* scene, * grouporlight;
		scene = doc.FirstChildElement();
		if (scene) {
			grouporlight = scene->FirstChildElement();
			while (grouporlight) {
				if (grouporlight->Name() == group) {
					ler_group(grouporlight);

				}
				if (grouporlight->Name() == lights) {
					;
					ler_lights(grouporlight);
				}
				grouporlight = grouporlight->NextSiblingElement();
			}
		}
	}
	else {
		break_xml = true;
		switch (flag) {
		case(1):
			printf("ERROR: XML_NO_ATTRIBUTE\n");
			break;
		case(2):
			printf("ERROR: XML_WRONG_ATTRIBUTE_TYPE\n");
			break;
		case(3):
			printf("ERROR: XML_ERROR_FILE_NOT_FOUND\n");
			break;
		case(4):
			printf("ERROR: XML_ERROR_FILE_COULD_NOT_BE_OPENED\n");
			break;
		case(5):
			printf("ERROR: XML_ERROR_FILE_READ_ERROR\n");
			break;
		case(6):
			printf("ERROR: XML_ERROR_PARSING_ELEMENT\n");
			break;
		case(7):
			printf("ERROR: XML_ERROR_PARSING_ATTRIBUTE\n");
			break;
		case(8):
			printf("ERROR: XML_ERROR_PARSING_TEXT\n");
			break;
		case(9):
			printf("ERROR: XML_ERROR_PARSING_CDATA\n");
			break;
		case(10):
			printf("ERROR: XML_ERROR_PARSING_COMMENT\n");
			break;
		case(11):
			printf("ERROR: XML_ERROR_PARSING_DECLARATION\n");
			break;
		case(12):
			printf("ERROR: XML_ERROR_PARSING_UNKNOWN\n");
			break;
		case(13):
			printf("ERROR: XML_ERROR_EMPTY_DOCUMENT\n");
			break;
		case(14):
			printf("ERROR: XML_ERROR_MISMATCHED_ELEMENT\n");
			break;
		case(15):
			printf("ERROR: XML_ERROR_PARSING\n");
			break;
		case(16):
			printf("ERROR: XML_CAN_NOT_CONVERT_TEXT\n");
			break;
		case(17):
			printf("ERROR: XML_NO_TEXT_NODE\n");
			break;
		case(18):
			printf("ERROR: XML_ELEMENT_DEPTH_EXCEEDED\n");
			break;
		default:
			printf("ERROR: UNKNOWN ERROR\n");
			break;
		}
	}
}

//Prepara as texturas a serem utilizadas
int carregaTextura(std::string s) {
	unsigned int t, tw, th;
	unsigned char* texData;
	unsigned int texID;

	ilGenImages(1, &t);
	ilBindImage(t);
	ilLoadImage((ILstring)s.c_str());
	tw = ilGetInteger(IL_IMAGE_WIDTH);
	th = ilGetInteger(IL_IMAGE_HEIGHT);
	ilConvertImage(IL_RGBA, IL_UNSIGNED_BYTE);
	texData = ilGetData();

	glGenTextures(1, &texID);

	glBindTexture(GL_TEXTURE_2D, texID);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, tw, th, 0, GL_RGBA, GL_UNSIGNED_BYTE, texData);
	glGenerateMipmap(GL_TEXTURE_2D);

	glBindTexture(GL_TEXTURE_2D, 0);

	return texID;

}

//Prepara o programa para utilizar funcionalidades de luz e texturas 
void initGL() {
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT);
	glFrontFace(GL_CW);

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
	
	glEnable(GL_LIGHTING);
	glEnable(GL_TEXTURE_2D);
	glewInit();
	lighinit();
	ilInit();
	ilEnable(IL_ORIGIN_SET);
	ilOriginFunc(IL_ORIGIN_LOWER_LEFT);
}

int main(int argc, char** argv) {
	xml = argv[1];
	leitor_xml(xml);
	if (break_xml) return 1;
	vector<string>argf;
	// put GLUT’s init here
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(1000, 1000);
	glutCreateWindow("Fase_4");

	// put callback registry here
	glutKeyboardFunc(input);
	glutMouseFunc(carregar_rato);
	glutMotionFunc(rato);
	glutReshapeFunc(changeSize);
	glutIdleFunc(renderScene);
	glutDisplayFunc(renderScene);
	timebase = glutGet(GLUT_ELAPSED_TIME);
	initGL();
	
	for (int i = 0; i < files.size(); i++) {
		if (files[i].find("teste") == 0) {
			argf = splitstring(files[i], ',', argf);
			if (argf[1] == "T") {
				argf[2] = to_string(carregaTextura(argf[2]));
				files[i] = "";
				for (int x = 0; x < argf.size() - 1; x++) {
					files[i] = files[i] + argf[x] + ",";
				}
				files[i] += argf[argf.size()-1];
			}
		}
		argf.clear();
	}
	glGenBuffers(10, vertexbuffers);
	glGenBuffers(10, normalbuffers);
	glGenBuffers(10, texturebuffers);
	desenha(files);
	for (int i = 0; i < files.size(); i++) {
		if (files[i].find("teste") == 0) {
			argf = splitstring(files[i], ',', argf);
			triangleCount += verticeCount[qual_objeto(argf[0])];
			argf.clear();
		}
	}
	triangleCount /= 3;
	printf("Botao_Esquerdo -> Mexer camara\nBotao_Direito -> Ajustar zoom (apenas para Third Person Camera)\nBotao_Meio -> Selecionar Objeto (apenas para sistemasolar.xml)\n");
	printf("f -> Desenhar Triangulos\nl -> Desenhar Linhas\np -> Desenhar Pontos\nn -> Desenhar Normais\nx -> Aumentar Velocidade\nz -> Diminuir Velocidade\n");
	printf("c -> Alterar Camara (First/Third Person Camera)\n");
	printf("w -> Mover para a Frente (FPC)\ns -> Mover para Tras (FPC)\na -> Mover para a Esquerda (FPC)\nd -> Mover para a Direita (FPC)\n");
	printf("e -> Subir (FPC)\nq -> Descer(FPC)\n");
	// enter GLUT’s main cycle
	glutMainLoop();
	return 1;
}