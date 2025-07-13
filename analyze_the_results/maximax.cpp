#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define R 1000								// �ЂƂ̐^�̏d�v�x�̑g�ɑ΂���v���O�����̌J��Ԃ���						
#define Generated_method 3				//������@ a3 or a4
#define N 6								// �]����̐�4-8�̒l�@���������͔z��̒�`�Ɏg�����߃��[�v���Ă��Ȃ�
#define Alt 5							// ��ֈĐ�
#define M 100							// ���p�s��̌�
#define bM 1
#define Rt 5							// �^�̏d�v�x��^����� A~E
#define Method_N 2						// ����@��@�̌�



int evaluation_num = N - 4;				//����unameLIst�̔z��ɓ���邽�߂̐���
char eval_num_filepath[5][31] = { "\\N=4" ,"\\N=5" ,"\\N=6" ,"\\N=7+" ,"\\N=8+" };//�]����̐��͂����������肷��z��
char util_num_filepath[5][31] = { "\\N=4" ,"\\N=5" ,"\\N=6_M=5" ,"\\N=7+" ,"\\N=8+" };//�]����̐��͂����������肷��z��


const double PI = 3.14159265358979;
double epsi = 0.000001;					// ���AHP�̃C�v�V����

char true_weight_list[Rt][8] = { "A","B","C","D","E" };		//A-E
//char method_name_list[31][64] = { "\\WMIN","\\WWMIN","\\MMRW","\\AMRW","\\E-MMRW","\\G-MMRW","\\DMIN","\\MMRD","\\AMRD","\\E-MMRD","\\G-MMRD","\\MMRLD","\\AMRLD",
//"\\EV","\\GM","\\E-WMIN","\\G-WMIN","\\E-DMIN","\\G-DMIN","\\E-AMRW","\\G-AMRW","\\E-AMRD","\\G-AMRD","\\MMRWW","\\AMRWW",
//"\\E-MMRWW","\\G-MMRWW","\\E-AMRWW","\\G-AMRWW","\\E-WWMIN","\\G-WWMIN" };//��@�̖��O
//
//char file_method[31][64] = {"\\MSW","\\MSWW","\\MMRW","\\AMRW","\\E-MMRW","\\G-MMRW","\\MSD","\\MMRD","\\AMRD","\\E-MMRD","\\G-MMRD",
//"","","\\EV","\\GM","\\E-MSW","\\G-MSW","\\E-MSD","\\G-MSD","\\E-AMRW","\\G-AMRW","\\E-AMRD","\\G-AMRD","\\MMRWW",
//"\\AMRWW","\\E-MMRWW","\\G-MMRWW","\\E-AMRWW","\\G-AMRWW","\\E-MSWW","\\G-MSWW"}; //�e��@�̃t�@�C���̖��O
//// �g�p���鐄��@��ID���X�g
//int method_id_list[Method_N] = { 0,1,2,3,4,5,6,7,8,9,10,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30 };


//char file_method[16][64] = { "E-AMRw", "E-MMRw", "G-AMRw", "G-MMRw", "eAMRw", "eAMRwc", "eMMRw", "eMMRwc", "gAMRw", "gAMRwc", "gMMRw", "gMMRwc", "lAMRw", "lAMRwc", "lMMRw", "lMMRwc" };

void BubSort_Ascend(double x[][2], int item);										// �\�[�g ����	Morm = 0;maximin�̂Ƃ�
void BubSort_Descend(double x[][2], int item);										//�\�[�g�@�~��  Morm = 1;Maximax�̂Ƃ�
void maximin(double u[Alt][N], double wR[N], double wL[N], double totalU[Alt], double z[Alt][N], int perm[Alt][N], int star[Alt]);	// �ŏ��������p�l�̌v�Z
int combi(int n, int m);														// �y�A�̔ԍ�
void arc_combi(int comb, int n[2]);

int main(void) {
	char max_or_min[2][64] = {"maximin","Maximax"};
	//maximin,Maximax���`����
	for (int Morm = 0;Morm < 2;Morm++) {
		//u1 or u2�̓��
		for (int utility = 1;utility < 3;utility++) {
			int utility_num = utility - 1;
			char which_utility[2][31] = { "\\u1","\\u2" };
			int f;
			int md, rt;
			//char method_name_list[16][64] = { "\\E-AMRw", "\\E-MMRw", "\\G-AMRw", "\\G-MMRw", "\\eAMRw", "\\eAMRwc", "\\eMMRw", "\\eMMRwc", "\\gAMRw", "\\gAMRwc", "\\gMMRw", "\\gMMRwc", "\\lAMRw", "\\lAMRwc", "\\lMMRw", "\\lMMRwc" };
			char method_name_list[2][64] = { "\\AMRwc","\\MMRwc" };
			char file_method[2][64] = { "AMRwc","MMRwc" };
			double trueW[2 * N], true_wR[N], true_wL[N];

			double ukari[Alt * M][8];


			FILE* fp_uM, * fp_iw;	// ���̓t�@�C��

			char filepath[256];

			// �^�̋�ԏd�v�xsetting
			for (rt = 0; rt < Rt; rt++) {

				printf("u%d:%s\n", utility, true_weight_list[rt]);

				// �^�̋�ԏd�v�x�ǂݍ���
				sprintf(filepath, "..\\data\\true_interval_weight_set%s\\%s\\Given_interval_weight.csv", eval_num_filepath[evaluation_num], true_weight_list[rt]);
				fp_iw = fopen(filepath, "r");
				if (fp_iw == nullptr) {
					printf("�G���[: �t�@�C�����J���܂���ł��� (fp_iw): %s\n", filepath);
				}
				printf("�J����");
				for (f = 0; f < 2 * N; f++) {
					fscanf(fp_iw, "%lf,", &trueW[f]);
				}
				for (f = 0; f < 2 * N; f++) {
					if (f % 2 == 0) true_wL[f / 2] = trueW[f];
					else true_wR[f / 2] = trueW[f];
				}
				fclose(fp_iw);

				// ���ӌ��p�lu�ǂݍ���
				sprintf(filepath, "..\\data\\���p�l�s��%s%s\\u.csv", which_utility[utility_num], util_num_filepath[evaluation_num]);
				fp_uM = fopen(filepath, "r");
				if (fp_uM == nullptr) {
					printf("�G���[: �t�@�C�����J���܂���ł��� (fp_iw): %s\n", filepath);
				}
				printf("�J����");
				for (f = 0; f < Alt * M; f++) {
					fscanf(fp_uM, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &ukari[f][0], &ukari[f][1], &ukari[f][2], &ukari[f][3], &ukari[f][4], &ukari[f][5], &ukari[f][6], &ukari[f][7]);
				}
				fclose(fp_uM);

#pragma omp parallel for            //���s�������邽�߂̕��͂��̂��߂ɂ��̒��ł����ϐ��͂��̒��Œ�`����
				for (md = 0; md < Method_N; md++) {	// �e����@�̃��[�v
					int thread_id = omp_get_thread_num();
					printf("Thread %d is processing iteration %d\n", thread_id, md);
					int i, j, l, n, m, rp;			//���[�v�p�ϐ�
					FILE* fp_cs, * fp_rc, * fp_tc;	// �o�̓t�@�C��
					FILE* fp_ew;
					int k[2];
					int method_id;					// ����@�̃��[�v�J�E���^
					double u[Alt][N];				// ���ӌ��p�lu[Alt][N]

					double wR[N], wL[N];				// �����ԏd�v�x
					double wkari[R][2 * N];			// �����ԏd�v�x
					double wsum;
					double yR[N], yL[N];			// Y = [t*wL, t*wR]
					double z[Alt][N];				// �������p�l�ŗp����d�v�x�@z in Y
					double yR2[N], yL2[N];			// Y' = [r*yL, r*yR]

					double totalU[Alt];				// �������p�l
					double totalU2[Alt];			// �������p�l'

					double t_snap;					// t�܂�Ȃ���_
					double t[101][3];				// t�ۊǗp�@�������܂�
					double tL, tR;
					double sL, sR;					// �^�̋�ԏd�v�x�̏ꍇ
					double s_snap;					// s�܂�Ȃ���_
					double sh[101][3];				// s�ۊǗp�@�������܂�
					double alphaR, alphaL;

					int perm[Alt][N];				// ���ӌ��p�l�����Ԗڂɑ傫���v�f��
					double udp[N][2];				// �����ֈĂ̊e��̌��p�ƌ��X�ǂ̊�����L�^
					int star[Alt];					// yi^L < zi < yi^R �ƂȂ�i�̒l
					int star2[Alt];
					double Sl, S[100];
					double r, r2[100][2];
					int it_t, it_s;
					int count1;
					int flag1;
					int flag2;						// �덷����
					double flag3[Alt];
					double swap1;
					int swap2;
					int rank1[100][Alt], rank2[100][Alt];
					double rank1c_s[100], rank2c_t[100];	// �����ω�s,t
					int count_rank1, count_rank2;
					double tbud[Alt][Alt], bud[Alt][Alt];
					int countC10[100][100];

					double store1, store2;
					char method_name[64];

					int id;
					char bff[1024];
					char* pbff, * ebff;
					char filepath_o[256];
					// ����@
					method_id =md;

					//// �����ԏd�v�x�̓ǂ݂�
					sprintf(filepath_o, "..\\data\\Simp\\Simp%s\\a%d\\%s%s\\Simp.csv", eval_num_filepath[evaluation_num], Generated_method, true_weight_list[rt], method_name_list[method_id]);
					fp_ew = fopen(filepath_o, "r");
					if (fp_ew == nullptr) {
						printf("�G���[: �t�@�C�����J���܂���ł��� (fp_iw): %s\n", filepath_o);
					}
					fgets(bff, 1024, fp_ew);
					fgets(bff, 1024, fp_ew);
					fgets(bff, 1024, fp_ew);
					// �����܂�Simp.csv�̃w�b�_�[����, �ȍ~, �����ԏd�v�x��ǂݍ���
					for (i = 0; i < R; i++) {
						fgets(bff, 1024, fp_ew);
						pbff = bff;
						id = strtod(pbff, &pbff);
						pbff++;
						for (j = 0; j < 2 * N; j++) {
							wkari[i][j] = strtod(pbff, &pbff);
							pbff++;
						}
						wsum = strtod(pbff, &pbff);
						pbff++;
					}
					fclose(fp_ew);


					// �o�̓t�@�C��
					sprintf(filepath_o, "..\\data\\a%d\\%s%s%s\\%s%s\\%s_count_pairs_in_squares_%d_u%d.csv", Generated_method, max_or_min[Morm], which_utility[utility_num], eval_num_filepath[evaluation_num], true_weight_list[rt], method_name_list[method_id], max_or_min[Morm],R,utility);
					fp_cs = fopen(filepath_o, "w");

					/*
					sprintf(filepath_o, "%s%s", method_name, "\\ranking_change.csv");
					fp_rc = fopen(filepath_o, "w");

					sprintf(filepath_o, "%s%s", method_name, "\\t_change.csv");
					fp_tc = fopen(filepath_o, "w");
					*/
					for (m = 0; m < M; m++) {

						// ���ӌ��p�l�̌Ăяo��
						for (i = 0; i < Alt; i++) {
							for (j = 0; j < N; j++) {
								u[i][j] = ukari[Alt * m + i][j];
							}
						}

						// ���ӌ��p�l����
						for (j = 0; j < Alt; j++) {
							for (i = 0; i < N; i++) {
								udp[i][0] = u[j][i];	// i�͂ǂ̊����\��
								udp[i][1] = i;
							}
							Morm == 0 ? BubSort_Ascend(udp, N) : BubSort_Descend(udp, N);
							for (i = 0; i < N; i++) {
								perm[j][i] = (int)udp[i][1];
							}
						}

						////�^�̋�ԏd�v�x�ɂ����Đ��K���𖞂����͈�sL,sR�̌v�Z
						alphaR = 10, alphaL = 0;	// ������
						for (i = 0; i < N; i++) {
							store1 = true_wL[i];
							store2 = true_wR[i];
							for (j = 0; j < N; j++) {
								if (j != i) {
									store1 += true_wR[j];
									store2 += true_wL[j];
								}
							}
							if (store1 < alphaR) alphaR = store1;
							if (store2 > alphaL) alphaL = store2;
						}
						sR = 1 / alphaL;
						sL = 1 / alphaR;
						if (sL > 1.0) sL = 1.0;
						if (sR < 1.0) sR = 1.0;
						// ������
						for (i = 0; i < 100; i++) {
							for (j = 0; j < 3; j++) {
								sh[i][j] = 0;
							}
						}
						for (i = 0; i < 100; i++) {
							rank1c_s[i] = 0;
						}
						it_s = 0;
						sh[it_s][0] = sR; sh[it_s][1] = 0; sh[it_s++][2] = 0;
						s_snap = sR;
						count_rank1 = 0;
						rank1c_s[count_rank1] = sR;
						for (i = 0; i < 100; i++) {
							for (j = 0; j < Alt; j++) {
								rank1[i][j] = 0;
							}
						}

						while (s_snap >= sL) {

							// t*w^L,t*w^R
							for (i = 0; i < N; i++) {
								yR[i] = true_wR[i] * s_snap;
								yL[i] = true_wL[i] * s_snap;
							}

							// �}�L�V�~���������p�l
							maximin(u, yR, yL, totalU, z, perm, star);

							if (it_s == 1) { // 1��ڂ���rank1���v�Z
								//printf("s:it_t=1\n");
								for (i = 0; i < Alt; i++) {
									flag3[i] = totalU[i];
									rank1[count_rank1][i] = i;
								}

								for (i = 0; i < Alt; ++i) {
									for (j = i + 1; j < Alt; ++j) {
										if (flag3[i] < flag3[j]) {
											swap1 = flag3[i];
											flag3[i] = flag3[j];
											flag3[j] = swap1;

											swap2 = rank1[count_rank1][i];
											rank1[count_rank1][i] = rank1[count_rank1][j];
											rank1[count_rank1][j] = swap2;
										}
									}
								}
								count_rank1++;
							}

							//printf("U	:%f, %f, %f, %f, %f\n", totalU[0], totalU[1], totalU[2], totalU[3], totalU[4]);

							// ���ɌX�����ς��_�𒲂ׂ�
							store1 = 0; Sl = 10;
							for (j = 0; j < Alt; j++) {
								store1 = yR[star[j]] - z[j][star[j]];
								if (Sl > store1) {
									Sl = store1;
									flag1 = j;
								}
							}
							r = 1 / (1 + Sl);
							flag2 = 0;

							// �����덷�X�L�b�v
							if (r == 1) {
								r = r - epsi;
								flag2 = -1;
							}

							// tL��菬�����Ȃ�Ȃ��悤��
							if (r * s_snap < sL) {
								r = sL / s_snap;
								// printf("end\n");
							}

							// �������p�lU'���v�Z
							for (i = 0; i < N; i++) {
								yR2[i] = yR[i] * r;
								yL2[i] = yL[i] * r;
							}

							maximin(u, yR2, yL2, totalU2, z, perm, star2);
							//printf("U2	:%f, %f, %f, %f, %f\n", totalU2[0], totalU2[1], totalU2[2], totalU2[3], totalU2[4]);

							// ����ւ��_�̗L�� �S�y�A����
							// ������
							for (i = 0; i < 100; i++) {
								S[i] = 0; r2[i][0] = 0; r2[i][1] = 0;
							}
							count1 = 0;

							for (i = 0; i < Alt - 1; i++) {
								for (j = i + 1; j < Alt; j++) {
									if ((totalU[i] - totalU[j]) * (totalU2[i] - totalU2[j]) <= 0) {
										//printf("change %d\n", combi(i, j));
										S[count1] = (-1) * (totalU[i] - totalU[j]) / (u[i][star[i]] - u[j][star[j]]);
										r2[count1][0] = 1 / (1 + S[count1]);
										r2[count1++][1] = combi(i, j);
									}
								}
							}

							// r'����
							if (count1 > 0) {
								BubSort_Ascend(r2, count1);
							}

							for (i = count1 - 1; i >= 0; i--) {
								sh[it_s][0] = s_snap * r2[i][0];
								sh[it_s][1] = 1;		// ����
								sh[it_s][2] = r2[i][1];

								rank1c_s[count_rank1] = sh[it_s][0];

								arc_combi((int)r2[i][1], k);
								//printf("s:arccombi\n");
								for (j = 0; j < Alt; j++) {
									if (rank1[count_rank1 - 1][j] == k[0]) {
										rank1[count_rank1][j] = k[1];
									}
									else if (rank1[count_rank1 - 1][j] == k[1]) {
										rank1[count_rank1][j] = k[0];
									}
									else {
										rank1[count_rank1][j] = rank1[count_rank1 - 1][j];
									}
								}

								it_s++;
								count_rank1++;
							}

							if (sh[it_s - 1][0] > s_snap * r) {
								sh[it_s][0] = s_snap * r;
								sh[it_s][1] = flag2;			// �܂�
								sh[it_s][2] = flag1;
								it_s++;
								s_snap = sh[it_s - 1][0];
							}
							/*else if (sh[it_s - 1][0] == s_snap * r) {
								s_snap -= epsi;
								if (s_snap < sL) {
									s_snap = sL;
								}
								sh[it_s][0] = s_snap;
								sh[it_s][1] = -1;			// �덷�X�L�b�v
								sh[it_s][2] = flag1;
								it_s++;
								printf("s:s=sr\n");
								flag2 = 1;	// �덷�X�L�b�v
							}*/
							else {
								//printf("s_error_skip\n");
								rank1c_s[count_rank1] = sL;
								break;
							}

							if (s_snap == sL) {
								rank1c_s[count_rank1] = sL;
								break;
							}

							if (it_s >= 100) {
								printf("%s:error_roop:it_t = %d in md = %d,rp = %d,m = %d\n", true_weight_list[rt], it_t, md, rp, m);
								break;
							}
						}
						/*
						// �o��
						fprintf(fp_rc, "%d,true, %d\n", m, count_rank1);
						for (i = 0; i < Alt; i++) {
							fprintf(fp_rc, " , ,");
							for (j = 0; j < count_rank1; j++) {
								fprintf(fp_rc, "%d,", rank1[count_rank1 - 1 - j][i]);
							}
							fprintf(fp_rc, "\n");
						}
						fprintf(fp_tc, "%d, true, %d\n", m, it_s);
						for (i = 0; i < 3; i++) {
							fprintf(fp_tc, " , ,");
							for (j = 0; j < it_s; j++) {
								fprintf(fp_tc, "%f,", sh[it_s - 1 - j][i]);
							}
							fprintf(fp_tc, "\n");
						}
						*/

						// �����ԏd�v�x�̃��[�v
						for (rp = 0; rp < R; rp++) {

							// �����ԏd�v�x�̌Ăяo��
							for (i = 0; i < N; i++) {
								wL[i] = wkari[rp][2 * i];
								wR[i] = wkari[rp][2 * i + 1];
							}

							// �����ԏd�v�x�ɂ����Đ��K���𖞂����͈�tL,tR�̌v�Z
							alphaR = 10, alphaL = 0;	// ������
							for (i = 0; i < N; i++) {
								store1 = wL[i];
								store2 = wR[i];
								for (j = 0; j < N; j++) {
									if (j != i) {
										store1 += wR[j];
										store2 += wL[j];
									}
								}
								if (store1 < alphaR) alphaR = store1;
								if (store2 > alphaL) alphaL = store2;
							}
							tR = 1 / alphaL;
							tL = 1 / alphaR;
							if (tL > 1.0) tL = 1.0;
							if (tR < 1.0) tR = 1.0;
							// ������
							for (i = 0; i < 100; i++) {
								for (j = 0; j < 3; j++) {
									t[i][j] = 0;
								}
							}
							for (i = 0; i < 100; i++) {
								rank2c_t[i] = 0;
							}
							it_t = 0;
							t[it_t][0] = tR; t[it_t][1] = 0; t[it_t++][2] = 0;
							t_snap = tR;
							count_rank2 = 0;
							rank2c_t[count_rank2] = tR;
							for (i = 0; i < 100; i++) {
								for (j = 0; j < Alt; j++) {
									rank2[i][j] = 0;
								}
							}

							while (t_snap >= tL) {

								// t*w^L,t*w^R
								for (i = 0; i < N; i++) {
									yR[i] = wR[i] * t_snap;
									yL[i] = wL[i] * t_snap;
								}

								// �}�L�V�~���������p�l
								maximin(u, yR, yL, totalU, z, perm, star);

								if (it_t == 1) { // 1��ڂ���rank2���v�Z
									for (i = 0; i < Alt; i++) {
										flag3[i] = totalU[i];
										rank2[count_rank2][i] = i;
									}

									for (i = 0; i < Alt; ++i) {
										for (j = i + 1; j < Alt; ++j) {
											if (flag3[i] < flag3[j]) {
												swap1 = flag3[i];
												flag3[i] = flag3[j];
												flag3[j] = swap1;


												swap2 = rank2[count_rank2][i];
												rank2[count_rank2][i] = rank2[count_rank2][j];
												rank2[count_rank2][j] = swap2;
											}
										}
									}
									count_rank2++;
								}

								//printf("U	:%f, %f, %f, %f, %f\n", totalU[0], totalU[1], totalU[2], totalU[3], totalU[4]);

								// ���ɌX�����ς��_�𒲂ׂ�
								store1 = 0; Sl = 10;
								for (j = 0; j < Alt; j++) {
									store1 = yR[star[j]] - z[j][star[j]];
									if (Sl > store1) {
										Sl = store1;
										flag1 = j;
									}
								}
								r = 1 / (1 + Sl);
								flag2 = 0;

								// �����덷�X�L�b�v
								if (r == 1) {
									r = r - epsi;
									flag2 = -1;
								}

								// tL��菬�����Ȃ�Ȃ��悤��
								if (r * t_snap < tL) {
									r = tL / t_snap;
									// printf("end\n");
								}

								// �������p�lU'���v�Z
								for (i = 0; i < N; i++) {
									yR2[i] = yR[i] * r;
									yL2[i] = yL[i] * r;
								}

								maximin(u, yR2, yL2, totalU2, z, perm, star2);
								//printf("U2	:%f, %f, %f, %f, %f\n", totalU2[0], totalU2[1], totalU2[2], totalU2[3], totalU2[4]);

								// ����ւ��_�̗L�� �S�y�A����
								// ������
								for (i = 0; i < 100; i++) {
									S[i] = 0; r2[i][0] = 0; r2[i][1] = 0;
								}
								count1 = 0;

								for (i = 0; i < Alt - 1; i++) {
									for (j = i + 1; j < Alt; j++) {
										if ((totalU[i] - totalU[j]) * (totalU2[i] - totalU2[j]) <= 0) {
											//printf("change %d\n", combi(i, j));
											S[count1] = (-1) * (totalU[i] - totalU[j]) / (u[i][star[i]] - u[j][star[j]]);
											r2[count1][0] = 1 / (1 + S[count1]);
											r2[count1++][1] = combi(i, j);
										}
									}
								}

								// r'����
								if (count1 > 0) {
									BubSort_Ascend(r2, count1);
								}

								for (i = count1 - 1; i >= 0; i--) {
									t[it_t][0] = t_snap * r2[i][0];
									t[it_t][1] = 1;		// ����
									t[it_t][2] = r2[i][1];

									arc_combi((int)r2[i][1], k);
									for (j = 0; j < Alt; j++) {
										if (rank2[count_rank2 - 1][j] == k[0]) {
											rank2[count_rank2][j] = k[1];
										}
										else if (rank2[count_rank2 - 1][j] == k[1]) {
											rank2[count_rank2][j] = k[0];
										}
										else {
											rank2[count_rank2][j] = rank2[count_rank2 - 1][j];
										}
									}
									rank2c_t[count_rank2] = t[it_t][0];
									count_rank2++;
									it_t++;
								}

								if (t[it_t - 1][0] > t_snap * r) {
									t[it_t][0] = t_snap * r;
									t[it_t][1] = flag2;			// �܂�
									t[it_t][2] = flag1;
									it_t++;
									t_snap = t[it_t - 1][0];
								}/*
								else if (t[it_t - 1][0] == t_snap * r) {
									t_snap -= epsi;
									if (t_snap < tL) {
										t_snap = tL;
									}
									t[it_t][0] = t_snap;
									t[it_t][1] = flag2;			// �덷�X�L�b�v
									t[it_t][2] = flag1;
									it_t++;
									//printf("t:t=tr\n");
									flag2 = 1;	// �덷�X�L�b�v
								}*/
								else {
									//printf("error_skip\n");
									rank2c_t[count_rank2] = tL;
									break;
								}

								if (t_snap <= tL) {
									rank2c_t[count_rank2] = tL;
									break;
								}

								if (it_t >= 100) {
									it_t = 100;
									printf("%s:error_roop:it_t = %d in md = %d,rp = %d,m = %d\n", true_weight_list[rt], it_t, md, rp, m);
									break;
								}
							}

							/*for (i = 0; i < it_t; i++) {
								for (j = 0; j < N; j++) {
									yL[j] = t[i][0] * wL[j];
									yR[j] = t[i][0] * wR[j];
								}
								maximin(u, yR, yL, totalU2, z, perm, star2);
								//printf("%d(%1f) t=%f	:%f, %f, %f, %f, %f\n", i, t[i][1], t[i][0], totalU2[0], totalU2[1], totalU2[2], totalU2[3], totalU2[4]);

							}*/

							/*
							// �o��
							fprintf(fp_rc, "%d,%d, %d\n", m, rp, count_rank2);
							for (i = 0; i < Alt; i++) {
								fprintf(fp_rc, " , ,");
								for (j = 0; j < count_rank2; j++) {
									fprintf(fp_rc, "%d,", rank2[count_rank2 - 1 - j][i]);
								}
								fprintf(fp_rc, "\n");
							}
							fprintf(fp_tc, "%d, %d, %d\n", m, rp, it_t);
							for (i = 0; i < 3; i++) {
								fprintf(fp_tc, " , ,");
								for (j = 0; j < it_t; j++) {
									fprintf(fp_tc, "%f,", t[it_t - 1 - j][i]);
								}
								fprintf(fp_tc, "\n");
							}
							*/
							for (i = 0; i < 100; i++) {
								for (j = 0; j < 100; j++) {
									countC10[i][j] = 0;
								}
							}

							// �召�֌W�̃J�E���g
							for (i = count_rank1 - 1; i >= 0; i--) {
								//�^
								for (j = 0; j < Alt; j++) {
									for (l = 0; l < Alt; l++) {
										tbud[j][l] = -1;
									}
								}
								for (j = 0; j < Alt; j++) {
									tbud[j][j] = 0;
								}
								for (j = 0; j < Alt; j++) {
									for (l = 0; l < Alt; l++) {
										if (tbud[rank1[i][j]][l] == -1) {
											tbud[rank1[i][j]][l] = 1;
										}
										if (tbud[l][rank1[i][j]] == -1) {
											tbud[l][rank1[i][j]] = 0;
										}
									}
								}

								// ����
								for (n = count_rank2 - 1; n >= 0; n--) {
									for (j = 0; j < Alt; j++) {
										for (l = 0; l < Alt; l++) {
											bud[j][l] = -1;
										}
									}
									for (j = 0; j < Alt; j++) {
										bud[j][j] = 0;
									}
									for (j = 0; j < Alt; j++) {
										for (l = 0; l < Alt; l++) {
											if (bud[rank2[n][j]][l] == -1) {
												bud[rank2[n][j]][l] = 1;
											}
											if (bud[l][rank2[n][j]] == -1) {
												bud[l][rank2[n][j]] = 0;
											}
										}
									}
									for (j = 0; j < Alt; j++) {
										for (l = 0; l < Alt; l++) {
											if (j < l) {
												if (bud[j][l] == tbud[j][l]) {
													countC10[count_rank1 - 1 - i][count_rank2 - 1 - n]++;
												}
											}
										}
									}
									//printf("%d,", countC10[count_rank1 - 1 - i][count_rank2 - 1 - n]);
								}
								//printf("\n");
							}
							//printf("\n");

							fprintf(fp_cs, "%d,%d,%d,%d\n", m, rp, count_rank1, count_rank2);
							l = count_rank2;
							fprintf(fp_cs, " , ,%f,", rank2c_t[l--]);
							for (j = 0; j <= count_rank1; j++) {
								fprintf(fp_cs, "%f,", rank1c_s[count_rank1 - j]);
							}
							fprintf(fp_cs, "\n");
							for (i = 0; i < count_rank2; i++) {
								fprintf(fp_cs, " , ,%f,", rank2c_t[l--]);
								for (j = 0; j < count_rank1; j++) {
									fprintf(fp_cs, "%d,", countC10[count_rank1 - 1 - j][count_rank2 - 1 - i]);
								}
								fprintf(fp_cs, "\n");
							}
						}
					}

					fclose(fp_cs);
					//fclose(fp_rc);
					//fclose(fp_tc);
				}

			}
		}
	}
	return 0;
}

/* maximin:�O�̗v�f�̕����傫����Ό�������*/
void BubSort_Ascend(double x[][2], int item) // �����\�[�g
{
	int i, j;
	double temp, temp1;

	for (i = 0; i < item - 1; i++) {
		for (j = item - 1; j > i; j--) {
			if( x[j - 1][0] > x[j][0])  {  
				

				temp = x[j][0];
				temp1 = x[j][1];

				x[j][0] = x[j - 1][0];
				x[j][1] = x[j - 1][1];

				x[j - 1][0] = temp;
				x[j - 1][1] = temp1;
			}
		}
	}
}
/* maximax:�O�̗v�f�̕�������������Ό�������*/
void BubSort_Descend(double x[][2], int item) // �~���\�[�g
{
	int i, j;
	double temp, temp1;

	for (i = 0; i < item - 1; i++) {
		for (j = item - 1; j > i; j--) {
			if (x[j - 1][0] < x[j][0]) {
				

				temp = x[j][0];
				temp1 = x[j][1];

				x[j][0] = x[j - 1][0];
				x[j][1] = x[j - 1][1];

				x[j - 1][0] = temp;
				x[j - 1][1] = temp1;
			}
		}
	}
}

void maximin(double u[Alt][N], double wR[N], double wL[N], double totalU[Alt], double z[Alt][N], int perm[Alt][N], int star[Alt])
{
	int i, j, k;
	int it;
	double cap;				// �ϐ��ł���d�v�x�̑��a��1�ƂȂ�悤�ɂ��邽��

	for (k = 0;k < Alt;k++) {

		//Process1
		it = 0;
		cap = 0.0;
		for (i = 0;i < N;i++) {
			cap = cap + wL[i];
		}

		//Process2
		while (cap + wR[perm[k][it]] - wL[perm[k][it]] <= 1) {
			z[k][perm[k][it]] = wR[perm[k][it]];
			cap = cap + wR[perm[k][it]] - wL[perm[k][it]];
			it = it + 1;
			if (it == N - 1) {
				break;
			}
		}

		//Process3
		z[k][perm[k][it]] = 1 - cap + wL[perm[k][it]];
		star[k] = perm[k][it];
		it = it + 1;

		//Process4
		while (it < N) {
			z[k][perm[k][it]] = wL[perm[k][it]];
			it = it + 1;
		}

		//�œK�l�̌v�Z
		totalU[k] = 0;
		for (i = 0;i < N;i++) {
			totalU[k] = totalU[k] + u[k][i] * z[k][i];
		}
	}
}

int combi(int n, int m)
{
	/*
	�Ⴆ�΁F��ֈ�4�Ńy�A�̔ԍ���U��Ƃ�
	0&1=0, 0&2=1, 0&3=2, 1&2=3, 1&3=4, 2&3=5�@�Ƃ����ӂ��ɐU���悤��
	(Alt) C 2�ʂ�
	*/

	int i, o, p;
	int comb = 0;		// �g�ݍ��킹�̔ԍ�

	if (n < m) {
		o = n + 1; p = m + 1;
	}
	else if (n > m) {
		o = m + 1; p = n + 1;
	}
	else {
		return -1;	// n,m���������ꍇ�̓G���[
	}

	for (i = 0; i < o - 1; i++) {
		comb = comb + (Alt - 1 - i);
	}

	comb = comb + p - o - 1;

	return comb;
}

void arc_combi(int comb, int n[2])
{
	int i = 0, j = comb + 1, o = 1, p = 2;

	while (j > 0) {
		if (j <= (Alt - o)) {
			p = j + o;
			break;
		}
		j = j - (Alt - o);
		o++;
	}

	n[0] = o - 1;
	n[1] = p - 1;
}