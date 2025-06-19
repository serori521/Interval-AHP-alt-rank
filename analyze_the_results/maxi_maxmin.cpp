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

#define R 1000								// �ЂƂ̐^�̏d�v�x�̑g�ɑ΂���v���O�����̌J��Ԃ���

#define Generated_method 3				//������@ a3 or a4
#define Alt 5							// ��ֈĐ�
#define M 100							// ���p�s��̌�
#define bM 1
#define Rt 5							// �^�̏d�v�x��^����� A~E
#define Method_N 16						// ����@��@�̌�


const double PI = 3.14159265358979;
double epsi = 0.000001;					// ���AHP�̃C�v�V����



int main(void) {

	// �^�̋�ԏd�v�xsetting
	omp_set_num_threads(4);
	char max_or_min[2][64] = { "maximin","Maximax" };
	//maximin,Maximax���`����
	#pragma omp parallel for collapse(2) schedule(static, 1)
	for (int Morm = 0; Morm < 2; Morm++) {
		for (int utility = 1;utility < 3;utility++){
			int utility_num = utility - 1;
			// �]����̐�4-8�̒l
			#pragma omp critical // printf�̏o�͂�������Ȃ��悤�ɕی�
			{
				// ���݂̑��X���b�h���A������ID�A�S�����Ă���g�ݍ��킹���o��
				printf("���X���b�h��: %d, �S���X���b�hID: %d  ==> [Morm=%d, utility=%d] ��������\n",
					omp_get_num_threads(), omp_get_thread_num(), Morm, utility);
			}

			for (int N = 6; N < 7;N++) {

				char which_utility[2][31] = { "\\u1","\\u2" };
				char true_weight_list[Rt][8] = { "A","B","C","D","E" };		//A-E
				char method_name_list[16][64] = { "\\E-AMRw", "\\E-MMRw", "\\G-AMRw", "\\G-MMRw", "\\eAMRw", "\\eAMRwc", "\\eMMRw", "\\eMMRwc", "\\gAMRw", "\\gAMRwc", "\\gMMRw", "\\gMMRwc", "\\lAMRw", "\\lAMRwc", "\\lMMRw", "\\lMMRwc" };
				char file_method[16][64] = { "E-AMRw", "E-MMRw", "G-AMRw", "G-MMRw", "eAMRw", "eAMRwc", "eMMRw", "eMMRwc", "gAMRw", "gAMRwc", "gMMRw", "gMMRwc", "lAMRw", "lAMRwc", "lMMRw", "lMMRwc" };
				//int method_id_list[Method_N] = {0};		// �g�p���鐄��@��ID���X�g
				//int method_id_list[Method_N] = { 0,1,2,3,4,5,6,7,8,9,10,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30 };
				int evaluation_num = N - 4;				//����unameLIst�̔z��ɓ���邽�߂̐���
				char given_weight[5][31] = { "\\N=4" ,"\\N=5" ,"\\N=6" ,"\\N=7+" ,"\\N=8+" };//�]����̐��͂����������肷��z��
				int i, j, k, l;
				int rt, md, m, rp;

				int count_rank1, count_rank2;
				double rank1c_s[100], rank2c_t[100];	// �����ω�����t,s�̒l
				int countC10[100][100];

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

						sprintf(filepath_o, "..\\data\\a3\\%s%s%s\\%s%s\\%s_count_pairs_in_squares_%d_u%d.csv", max_or_min[Morm], which_utility[utility_num],given_weight[evaluation_num], true_weight_list[rt], method_name_list[method_id], max_or_min[Morm], R, utility);
						printf("%s\n", filepath_o);
						if (filepath_o == NULL) {
							// �t�@�C�����J���Ȃ������ꍇ�̏���
							printf("�G���[: �t�@�C�����J���܂���ł����B�p�X���m�F���Ă��������B\n");
							printf("PATH: %s\n", filepath_o);
							continue; // md���[�v�̎��̌J��Ԃ��֐i�ށi�܂��� return �Ȃǂŏ����𒆒f�j
						}
						fp_cs = fopen(filepath_o, "r");
						// ����@���̏�����
						sum_rec = 0;
						sum_pre = 0;
						// ���ӌ��p�l���Ƃ̃��[�v
						for (m = 0; m < M; m++) {

							// �����ԏd�v�x�̃��[�v
							for (rp = 0; rp < R; rp++) {

								// �����ԏd�v�x���̏�����
								for (i = 0; i < 100; i++) {
									rank1c_s[i] = 0;
									rank2c_t[i] = 0;
									for (j = 0; j < 100; j++) {
										countC10[i][j] = 0;
									}
									flag1[i] = 0;
								}
								rec = 0;
								pre = 0;
								fscanf(fp_cs, "%d,%d,%d,%d\n", &k, &l, &count_rank1, &count_rank2);
								fscanf(fp_cs, " , ,%lf,", &rank2c_t[0]);
								for (j = 0; j <= count_rank1; j++) {
									fscanf(fp_cs, "%lf,", &rank1c_s[j]);
								}
								for (i = 0; i < count_rank2; i++) {
									fscanf(fp_cs, " , ,%lf,", &rank2c_t[i + 1]);
									for (j = 0; j < count_rank1; j++) {
										fscanf(fp_cs, "%d,", &countC10[i][j]);
									}
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
									
								}
								rec = rec / (double)count_rank1;

								sum_rec = sum_rec + rec;



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
								sum_pre = sum_pre + pre;
							}
						}

						sum_recall[rt][md] = sum_rec;
						sum_precision[rt][md] = sum_pre;

				
						fclose(fp_cs);

					}

				}
				sprintf(filepath_o, "..\\data\\%s_results%s%s�Č�������%d_u%d_%s.csv", max_or_min[Morm], which_utility[utility_num], given_weight[evaluation_num], R, utility, max_or_min[Morm]);
				fp_sum_recall = fopen(filepath_o, "w");
				fprintf(fp_sum_recall, "%s\n,A,B,C,D,E\n", filepath_o);

				sprintf(filepath_o, "..\\data\\%s_results%s%s�ŗǐ��𐔕���%d_u%d_%s.csv", max_or_min[Morm], which_utility[utility_num], given_weight[evaluation_num], R, utility, max_or_min[Morm]);
				fp_sum_precision = fopen(filepath_o, "w");
				fprintf(fp_sum_precision, "%s\n,A,B,C,D,E\n", filepath_o);

				//1�s�ڂ��L��
				//�f�[�^�̓���
				for (i = 0; i < Method_N; i++) {
					method_id = i;
					fprintf(fp_sum_precision, "%s,", file_method[method_id]);
					fprintf(fp_sum_recall, "%s,", file_method[method_id]);
					for (j = 0; j < Rt; j++) {
						fprintf(fp_sum_precision, "%f,", sum_precision[j][i]);
						fprintf(fp_sum_recall, "%f,", sum_recall[j][i]);
					}
					fprintf(fp_sum_precision, "\n");
					fprintf(fp_sum_recall, "\n");
				}
				fclose(fp_sum_recall);
				fclose(fp_sum_precision);

			}
		}
	}
	return 0;
}

