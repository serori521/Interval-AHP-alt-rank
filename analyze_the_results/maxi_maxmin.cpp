//********************************************************
// 区間重要度最適解集合全体でのペア一致数の解析
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

#define R 1000								// ひとつの真の重要度の組に対するプログラムの繰り返し回数

#define Generated_method 3				//生成手法 a3 or a4
#define Alt 5							// 代替案数
#define M 100							// 効用行列の個数
#define bM 1
#define Rt 5							// 真の重要度を与える回数 A~E
#define Method_N 16						// 推定法手法の個数


const double PI = 3.14159265358979;
double epsi = 0.000001;					// 区間AHPのイプシロン



int main(void) {

	// 真の区間重要度setting
	omp_set_num_threads(4);
	char max_or_min[2][64] = { "maximin","Maximax" };
	//maximin,Maximaxを定義する
	#pragma omp parallel for collapse(2) schedule(static, 1)
	for (int Morm = 0; Morm < 2; Morm++) {
		for (int utility = 1;utility < 3;utility++){
			int utility_num = utility - 1;
			// 評価基準の数4-8の値
			#pragma omp critical // printfの出力が混ざらないように保護
			{
				// 現在の総スレッド数、自分のID、担当している組み合わせを出力
				printf("総スレッド数: %d, 担当スレッドID: %d  ==> [Morm=%d, utility=%d] を処理中\n",
					omp_get_num_threads(), omp_get_thread_num(), Morm, utility);
			}

			for (int N = 6; N < 7;N++) {

				char which_utility[2][31] = { "\\u1","\\u2" };
				char true_weight_list[Rt][8] = { "A","B","C","D","E" };		//A-E
				char method_name_list[16][64] = { "\\E-AMRw", "\\E-MMRw", "\\G-AMRw", "\\G-MMRw", "\\eAMRw", "\\eAMRwc", "\\eMMRw", "\\eMMRwc", "\\gAMRw", "\\gAMRwc", "\\gMMRw", "\\gMMRwc", "\\lAMRw", "\\lAMRwc", "\\lMMRw", "\\lMMRwc" };
				char file_method[16][64] = { "E-AMRw", "E-MMRw", "G-AMRw", "G-MMRw", "eAMRw", "eAMRwc", "eMMRw", "eMMRwc", "gAMRw", "gAMRwc", "gMMRw", "gMMRwc", "lAMRw", "lAMRwc", "lMMRw", "lMMRwc" };
				//int method_id_list[Method_N] = {0};		// 使用する推定法のIDリスト
				//int method_id_list[Method_N] = { 0,1,2,3,4,5,6,7,8,9,10,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30 };
				int evaluation_num = N - 4;				//下のunameLIstの配列に入れるための数字
				char given_weight[5][31] = { "\\N=4" ,"\\N=5" ,"\\N=6" ,"\\N=7+" ,"\\N=8+" };//評価基準の数はいくつかを決定する配列
				int i, j, k, l;
				int rt, md, m, rp;

				int count_rank1, count_rank2;
				double rank1c_s[100], rank2c_t[100];	// 順序変化するt,sの値
				int countC10[100][100];

				double flag1[100];
				double rec;						// 再現率平均
				double rec_log;					// 対数重みづけ　再現率平均
				double pre;						// 正解率平均
				double pre_log;					// 対数重みづけ　正解率平均

				double sum_rec, sum_rec_log;	// 手法ごとの再現率平均の総和
				double sum_recall[Rt][Method_N], sum_recall_log[Rt][Method_N];
				double sum_pre, sum_pre_log;	// 手法ごとの正解率平均の総和
				double sum_precision[Rt][Method_N], sum_precision_log[Rt][Method_N];

				FILE* fp_cs;					// 入力ファイル
				FILE* fp_r, * fp_rl, * fp_mr, * fp_mrl;	// 再現率出力ファイル
				FILE* fp_p, * fp_pl, * fp_mp, * fp_mpl;	// 正解率出力ファイル
				FILE* fp_sum_recall, * fp_logsum_recall, * fp_sum_precision, * fp_logsum_precision;

				char filepath[256];
				char filepath_o[256];
				int method_id;					// 推定法のループカウンタ
				char method_name[64];

				for (rt = 0; rt < Rt; rt++) {

					printf("%d:%s\n", N, true_weight_list[rt]);

					// 推定法ごとのループ
					for (md = 0; md < Method_N; md++) {
						method_id = md;

						sprintf(filepath_o, "..\\data\\a3\\%s%s%s\\%s%s\\%s_count_pairs_in_squares_%d_u%d.csv", max_or_min[Morm], which_utility[utility_num],given_weight[evaluation_num], true_weight_list[rt], method_name_list[method_id], max_or_min[Morm], R, utility);
						printf("%s\n", filepath_o);
						if (filepath_o == NULL) {
							// ファイルが開けなかった場合の処理
							printf("エラー: ファイルが開けませんでした。パスを確認してください。\n");
							printf("PATH: %s\n", filepath_o);
							continue; // mdループの次の繰り返しへ進む（または return などで処理を中断）
						}
						fp_cs = fopen(filepath_o, "r");
						// 推定法毎の初期化
						sum_rec = 0;
						sum_pre = 0;
						// 周辺効用値ごとのループ
						for (m = 0; m < M; m++) {

							// 推定区間重要度のループ
							for (rp = 0; rp < R; rp++) {

								// 推定区間重要度毎の初期化
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
								////////////////// 再現率recall ////////////////////

								// 真における最大のマス
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



								////////////////// 正解率precision ////////////////////


								// 推定区間重要度毎の初期化
								for (i = 0; i < 100; i++) {
									flag1[i] = 0;
								}

								// 推定における最大のマス
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
				sprintf(filepath_o, "..\\data\\%s_results%s%s再現数平均%d_u%d_%s.csv", max_or_min[Morm], which_utility[utility_num], given_weight[evaluation_num], R, utility, max_or_min[Morm]);
				fp_sum_recall = fopen(filepath_o, "w");
				fprintf(fp_sum_recall, "%s\n,A,B,C,D,E\n", filepath_o);

				sprintf(filepath_o, "..\\data\\%s_results%s%s最良正解数平均%d_u%d_%s.csv", max_or_min[Morm], which_utility[utility_num], given_weight[evaluation_num], R, utility, max_or_min[Morm]);
				fp_sum_precision = fopen(filepath_o, "w");
				fprintf(fp_sum_precision, "%s\n,A,B,C,D,E\n", filepath_o);

				//1行目を記入
				//データの入力
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

