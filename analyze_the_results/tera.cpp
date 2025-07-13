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

#define R 10								// ひとつの真の重要度の組に対するプログラムの繰り返し回数
#define Generated_method 3				//生成手法 a3 or a4
#define Alt 5							// 代替案数
#define M 100							// 効用行列の個数
#define bM 1
#define Rt 5							// 真の重要度を与える回数 A~E
#define Method_N 2						// 推定法手法の個数


const double PI = 3.14159265358979;
double epsi = 0.000001;					// 区間AHPのイプシロン



int main(void) {

	// 真の区間重要度setting


	//u1 or u2という二択
	for (int utility = 2; utility < 3; utility++) {
		int utility_num = utility - 1;
		// 評価基準の数4-8の値
#pragma omp parallel for
		for (int N = 6; N < 7; N++) {

			char which_utility[2][31] = { "\\u1","\\u2" };
			char true_weight_list[Rt][8] = { "A","B","C","D","E" };		//A-E

			/*char method_name_list[16][64] = { "\\E-AMRw", "\\E-MMRw", "\\G-AMRw", "\\G-MMRw", "\\eAMRw", "\\eAMRwc", "\\eMMRw", "\\eMMRwc", "\\gAMRw", "\\gAMRwc", "\\gMMRw", "\\gMMRwc", "\\lAMRw", "\\lAMRwc", "\\lMMRw", "\\lMMRwc" };
			char file_method[16][64] = { "E-AMRw", "E-MMRw", "G-AMRw", "G-MMRw", "eAMRw", "eAMRwc", "eMMRw", "eMMRwc", "gAMRw", "gAMRwc", "gMMRw", "gMMRwc", "lAMRw", "lAMRwc", "lMMRw", "lMMRwc" };*/
			char method_name_list[2][64] = { "\\AMRwc","\\MMRwc" };
			char file_method[2][64] = { "AMRwc","MMRwc" };
			int evaluation_num = N - 4;				//下のunameLIstの配列に入れるための数字
			char given_weight[5][31] = { "\\N=4" ,"\\N=5" ,"\\N=6" ,"\\N=7+" ,"\\N=8+" };//評価基準の数はいくつかを決定する配列
			int i, j, k, l;
			int rt, md, m, rp;

			int count_rank1, count_rank2;
			double rank1c_s[100], rank2c_t[100];	// 順序変化するt,sの値
			double countC10[100][100];

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

					sprintf(filepath_o, "..\\data\\a3\\regret%s%s\\%s%s\\tera_minimax_regret_%d.csv", which_utility[utility_num], given_weight[evaluation_num], true_weight_list[rt], method_name_list[method_id], R);
					fp_cs = fopen(filepath_o, "r");



					// 推定法毎の初期化
					sum_rec = 0;
					sum_rec_log = 0;
					sum_pre = 0;
					sum_pre_log = 0;

					// 周辺効用値ごとのループ
					for (m = 0; m < M; m++) {

						// 推定区間重要度のループ
						for (rp = 0; rp < R; rp++) {

							// 推定区間重要度毎の初期化
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

							// 2行目: rank1c_s を読み込む
							// 先頭のカンマを読み込み、最初の double を読み込む (先頭のスペースを削除)
							fscanf(fp_cs, ",%lf", &rank1c_s[0]);

							// count_rank1 > 1 の場合、残りの double を読み込む
							for (j = 0; j < count_rank1 - 1; j++) {
								fscanf(fp_cs, ",%lf", &rank1c_s[j + 1]); // カンマに続けて double を読み込む
							}
							// 2行目の最後にある改行文字を読み飛ばす
							fscanf(fp_cs, "\n");

							// 3行目以降: rank2c_t と countC10 を読み込む (count_rank2 回繰り返す)
							for (i = 0; i < count_rank2; i++) {
								// 最初の double (rank2c_t) を読み込む (先頭のスペースを削除)
								fscanf(fp_cs, "%lf", &rank2c_t[i]);

								// countC10 の値を count_rank1 個読み込む
								for (j = 0; j < count_rank1; j++) {
									// カンマに続けて double (%lf) を読み込む ( %d から %lf に変更)
									fscanf(fp_cs, ",%lf", &countC10[i][j]);
								}
								// 各行の最後にある改行文字を読み飛ばす
								fscanf(fp_cs, "\n");
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
								//rec_log = rec_log + flag1[i] * log(rank1c_s[i + 1] / rank1c_s[i]);
							}
							rec = rec / (double)count_rank1;
							//rec_log = rec_log / log(rank1c_s[count_rank1] / rank1c_s[0]);


							sum_rec = sum_rec + rec;
							//sum_rec_log = sum_rec_log + rec_log;


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
			sprintf(filepath_o, "..\\data\\regret_results%s%stera再現数平均wc%d_u%d.csv", which_utility[utility_num], given_weight[evaluation_num], R, utility);
			fp_sum_recall = fopen(filepath_o, "w");
			fprintf(fp_sum_recall, "%s\n,A,B,C,D,E\n", filepath_o);

			//sprintf(filepath_o, "regret%s%s%s再現数加重平均%d_u%d.csv",  which_utility[utility_num], given_weight[evaluation_num], given_weight[evaluation_num], R,utility);
			//fp_logsum_recall = fopen(filepath_o, "w");
			//fprintf(fp_logsum_recall, "%s\n,A,B,C,D,E\n", filepath_o);

			sprintf(filepath_o, "..\\data\\regret_results%s%stera最良正解数平均wc%d_u%d.csv", which_utility[utility_num], given_weight[evaluation_num], R, utility);
			fp_sum_precision = fopen(filepath_o, "w");
			fprintf(fp_sum_precision, "%s\n,A,B,C,D,E\n", filepath_o);

			//sprintf(filepath_o, "regret%s%s%s最良正解数加重平均%d_u%d.csv",  which_utility[utility_num], given_weight[evaluation_num], given_weight[evaluation_num], R,utility);
			//fp_logsum_precision = fopen(filepath_o, "w");
			//fprintf(fp_logsum_precision, "%s\n,A,B,C,D,E\n", filepath_o);
			//1行目を記入
			//データの入力
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

