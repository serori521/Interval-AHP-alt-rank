//********************************************************
// ��ԏd�v�x�œK���W���S�̂ł̃y�A��v���̉��
// Author : Hayashi Akiko
// Update : 2023/1/25
//********************************************************

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define R 10								// �ЂƂ̐^�̏d�v�x�̑g�ɑ΂���v���O�����̌J��Ԃ���
#define Generated_method 3				//������@ a3 or a4
#define Alt 5							// ��ֈĐ�
#define M 100							// ���p�s��̌�
#define bM 1
#define Rt 5							// �^�̏d�v�x��^����� A~E
#define Method_N 2						// ����@��@�̌�


const double PI = 3.14159265358979;
double epsi = 0.000001;					// ���AHP�̃C�v�V����



int main(void) {

	// �^�̋�ԏd�v�xsetting


	//u1 or u2�Ƃ������
	for (int utility = 2; utility < 3; utility++) {
		int utility_num = utility - 1;
		// �]����̐�4-8�̒l
#pragma omp parallel for
		for (int N = 6; N < 7; N++) {

			char which_utility[2][31] = { "\\u1","\\u2" };
			char true_weight_list[Rt][8] = { "A","B","C","D","E" };		//A-E

			/*char method_name_list[16][64] = { "\\E-AMRw", "\\E-MMRw", "\\G-AMRw", "\\G-MMRw", "\\eAMRw", "\\eAMRwc", "\\eMMRw", "\\eMMRwc", "\\gAMRw", "\\gAMRwc", "\\gMMRw", "\\gMMRwc", "\\lAMRw", "\\lAMRwc", "\\lMMRw", "\\lMMRwc" };
			char file_method[16][64] = { "E-AMRw", "E-MMRw", "G-AMRw", "G-MMRw", "eAMRw", "eAMRwc", "eMMRw", "eMMRwc", "gAMRw", "gAMRwc", "gMMRw", "gMMRwc", "lAMRw", "lAMRwc", "lMMRw", "lMMRwc" };*/
			char method_name_list[2][64] = { "\\AMRwc","\\MMRwc" };
			char file_method[2][64] = { "AMRwc","MMRwc" };
			int evaluation_num = N - 4;				//����unameLIst�̔z��ɓ���邽�߂̐���
			char given_weight[5][31] = { "\\N=4" ,"\\N=5" ,"\\N=6" ,"\\N=7+" ,"\\N=8+" };//�]����̐��͂����������肷��z��
			int i, j, k, l;
			int rt, md, m, rp;

			int count_rank1, count_rank2;
			double rank1c_s[100], rank2c_t[100];	// �����ω�����t,s�̒l
			double countC10[100][100];

			double flag1[100];
			double rec;						// �Č�������
			double rec_log;					// �ΐ��d�݂Â��@�Č�������
			double pre;						// ���𗦕���
			double pre_log;					// �ΐ��d�݂Â��@���𗦕���

			double sum_rec, sum_rec_log;	// ��@���Ƃ̍Č������ς̑��a
			double sum_recall[Rt][Method_N], sum_recall_log[Rt][Method_N];
			double sum_pre, sum_pre_log;	// ��@���Ƃ̐��𗦕��ς̑��a
			double sum_precision[Rt][Method_N], sum_precision_log[Rt][Method_N];

			FILE* fp_cs;					// ���̓t�@�C��
			FILE* fp_r, * fp_rl, * fp_mr, * fp_mrl;	// �Č����o�̓t�@�C��
			FILE* fp_p, * fp_pl, * fp_mp, * fp_mpl;	// ���𗦏o�̓t�@�C��
			FILE* fp_sum_recall, * fp_logsum_recall, * fp_sum_precision, * fp_logsum_precision;

			char filepath[256];
			char filepath_o[256];
			int method_id;					// ����@�̃��[�v�J�E���^
			char method_name[64];
			for (rt = 0; rt < Rt; rt++) {

				printf("%d:%s\n", N, true_weight_list[rt]);

				// ����@���Ƃ̃��[�v
				for (md = 0; md < Method_N; md++) {
					method_id = md;

					sprintf(filepath_o, "..\\data\\a3\\regret%s%s\\%s%s\\tera_minimax_regret_%d.csv", which_utility[utility_num], given_weight[evaluation_num], true_weight_list[rt], method_name_list[method_id], R);
					fp_cs = fopen(filepath_o, "r");



					// ����@���̏�����
					sum_rec = 0;
					sum_rec_log = 0;
					sum_pre = 0;
					sum_pre_log = 0;

					// ���ӌ��p�l���Ƃ̃��[�v
					for (m = 0; m < M; m++) {

						// �����ԏd�v�x�̃��[�v
						for (rp = 0; rp < R; rp++) {

							// �����ԏd�v�x���̏�����
							for (i = 0; i < 100; i++) {
								rank1c_s[i] = 0;
								rank2c_t[i] = 0;
								for (j = 0; j < 100; j++) {
									countC10[i][j] = 0.0;
								}
								flag1[i] = 0;
							}
							rec = 0;
							rec_log = 0;
							pre = 0;
							pre_log = 0;

							fscanf(fp_cs, "%d,%d,%d,%d\n", &k, &l, &count_rank1, &count_rank2);

							// 2�s��: rank1c_s ��ǂݍ���
							// �擪�̃J���}��ǂݍ��݁A�ŏ��� double ��ǂݍ��� (�擪�̃X�y�[�X���폜)
							fscanf(fp_cs, ",%lf", &rank1c_s[0]);

							// count_rank1 > 1 �̏ꍇ�A�c��� double ��ǂݍ���
							for (j = 0; j < count_rank1 - 1; j++) {
								fscanf(fp_cs, ",%lf", &rank1c_s[j + 1]); // �J���}�ɑ����� double ��ǂݍ���
							}
							// 2�s�ڂ̍Ō�ɂ�����s������ǂݔ�΂�
							fscanf(fp_cs, "\n");

							// 3�s�ڈȍ~: rank2c_t �� countC10 ��ǂݍ��� (count_rank2 ��J��Ԃ�)
							for (i = 0; i < count_rank2; i++) {
								// �ŏ��� double (rank2c_t) ��ǂݍ��� (�擪�̃X�y�[�X���폜)
								fscanf(fp_cs, "%lf", &rank2c_t[i]);

								// countC10 �̒l�� count_rank1 �ǂݍ���
								for (j = 0; j < count_rank1; j++) {
									// �J���}�ɑ����� double (%lf) ��ǂݍ��� ( %d ���� %lf �ɕύX)
									fscanf(fp_cs, ",%lf", &countC10[i][j]);
								}
								// �e�s�̍Ō�ɂ�����s������ǂݔ�΂�
								fscanf(fp_cs, "\n");
							}

							////////////////// �Č���recall ////////////////////

							// �^�ɂ�����ő�̃}�X
							for (i = 0; i < count_rank1; i++) {
								flag1[i] = countC10[0][i];
								for (j = 1; j < count_rank2; j++) {
									if (countC10[j][i] > flag1[i]) {
										flag1[i] = countC10[j][i];
									}
								}
							}

							for (i = 0; i < count_rank1; i++) {
								rec = rec + flag1[i];
								//rec_log = rec_log + flag1[i] * log(rank1c_s[i + 1] / rank1c_s[i]);
							}
							rec = rec / (double)count_rank1;
							//rec_log = rec_log / log(rank1c_s[count_rank1] / rank1c_s[0]);


							sum_rec = sum_rec + rec;
							//sum_rec_log = sum_rec_log + rec_log;


							////////////////// ����precision ////////////////////

							// �����ԏd�v�x���̏�����
							for (i = 0; i < 100; i++) {
								flag1[i] = 0;
							}

							// ����ɂ�����ő�̃}�X
							for (i = 0; i < count_rank2; i++) {
								flag1[i] = countC10[i][0];
								for (j = 1; j < count_rank1; j++) {
									if (countC10[i][j] > flag1[i]) {
										flag1[i] = countC10[i][j];
									}
								}
							}

							for (i = 0; i < count_rank2; i++) {
								pre = pre + flag1[i];
							}
							pre = pre / (double)count_rank2;

							//// precision log
							//if ((count_rank2 == 1) || (rank2c_t[count_rank2] - rank2c_t[0] < epsi)) {
							//	pre_log = flag1[0];
							//}
							//else {
							//	for (i = 0; i < count_rank2; i++) {
							//		pre_log = pre_log + flag1[i] * log(rank2c_t[i + 1] / rank2c_t[i]);
							//	}
							//	pre_log = pre_log / log(rank2c_t[count_rank2] / rank2c_t[0]);
							//}


							sum_pre = sum_pre + pre;
							//sum_pre_log = sum_pre_log + pre_log;
						}
					}

					sum_recall[rt][md] = sum_rec;
					//sum_recall_log[rt][md] = sum_rec_log;
					sum_precision[rt][md] = sum_pre;
					//sum_precision_log[rt][md] = sum_pre_log;


					fclose(fp_cs);

				}

			}
			sprintf(filepath_o, "..\\data\\regret_results%s%stera�Č�������wc%d_u%d.csv", which_utility[utility_num], given_weight[evaluation_num], R, utility);
			fp_sum_recall = fopen(filepath_o, "w");
			fprintf(fp_sum_recall, "%s\n,A,B,C,D,E\n", filepath_o);

			//sprintf(filepath_o, "regret%s%s%s�Č������d����%d_u%d.csv",  which_utility[utility_num], given_weight[evaluation_num], given_weight[evaluation_num], R,utility);
			//fp_logsum_recall = fopen(filepath_o, "w");
			//fprintf(fp_logsum_recall, "%s\n,A,B,C,D,E\n", filepath_o);

			sprintf(filepath_o, "..\\data\\regret_results%s%stera�ŗǐ��𐔕���wc%d_u%d.csv", which_utility[utility_num], given_weight[evaluation_num], R, utility);
			fp_sum_precision = fopen(filepath_o, "w");
			fprintf(fp_sum_precision, "%s\n,A,B,C,D,E\n", filepath_o);

			//sprintf(filepath_o, "regret%s%s%s�ŗǐ��𐔉��d����%d_u%d.csv",  which_utility[utility_num], given_weight[evaluation_num], given_weight[evaluation_num], R,utility);
			//fp_logsum_precision = fopen(filepath_o, "w");
			//fprintf(fp_logsum_precision, "%s\n,A,B,C,D,E\n", filepath_o);
			//1�s�ڂ��L��
			//�f�[�^�̓���
			for (i = 0; i < Method_N; i++) {
				method_id = i;
				fprintf(fp_sum_precision, "%s,", file_method[method_id]);
				//fprintf(fp_logsum_precision, "%s,", file_method[method_id]);
				fprintf(fp_sum_recall, "%s,", file_method[method_id]);
				//fprintf(fp_logsum_recall, "%s,", file_method[method_id]);
				for (j = 0; j < Rt; j++) {
					fprintf(fp_sum_precision, "%f,", sum_precision[j][i]);
					//fprintf(fp_logsum_precision, "%f,", sum_precision_log[j][i]);
					fprintf(fp_sum_recall, "%f,", sum_recall[j][i]);
					//fprintf(fp_logsum_recall, "%f,", sum_recall_log[j][i]);
				}
				fprintf(fp_sum_precision, "\n");
				//fprintf(fp_logsum_precision, "\n");
				fprintf(fp_sum_recall, "\n");
				//fprintf(fp_logsum_recall, "\n");
			}
			fclose(fp_sum_recall);
			//fclose(fp_logsum_recall);
			fclose(fp_sum_precision);
			//fclose(fp_logsum_precision);
		}
	}
	return 0;
}

