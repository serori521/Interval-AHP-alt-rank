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
#define N 6								// 評価基準の数4-8の値　ここだけは配列の定義に使うためループしていない
#define Alt 5							// 代替案数
#define M 100							// 効用行列の個数
#define bM 1
#define Rt 5							// 真の重要度を与える回数 A~E
#define Method_N 2						// 推定法手法の個数



int evaluation_num = N - 4;				//下のunameLIstの配列に入れるための数字
char eval_num_filepath[5][31] = { "\\N=4" ,"\\N=5" ,"\\N=6" ,"\\N=7+" ,"\\N=8+" };//評価基準の数はいくつかを決定する配列
char util_num_filepath[5][31] = { "\\N=4" ,"\\N=5" ,"\\N=6_M=5" ,"\\N=7+" ,"\\N=8+" };//評価基準の数はいくつかを決定する配列


const double PI = 3.14159265358979;
double epsi = 0.000001;					// 区間AHPのイプシロン

char true_weight_list[Rt][8] = { "A","B","C","D","E" };		//A-E
//char method_name_list[31][64] = { "\\WMIN","\\WWMIN","\\MMRW","\\AMRW","\\E-MMRW","\\G-MMRW","\\DMIN","\\MMRD","\\AMRD","\\E-MMRD","\\G-MMRD","\\MMRLD","\\AMRLD",
//"\\EV","\\GM","\\E-WMIN","\\G-WMIN","\\E-DMIN","\\G-DMIN","\\E-AMRW","\\G-AMRW","\\E-AMRD","\\G-AMRD","\\MMRWW","\\AMRWW",
//"\\E-MMRWW","\\G-MMRWW","\\E-AMRWW","\\G-AMRWW","\\E-WWMIN","\\G-WWMIN" };//手法の名前
//
//char file_method[31][64] = {"\\MSW","\\MSWW","\\MMRW","\\AMRW","\\E-MMRW","\\G-MMRW","\\MSD","\\MMRD","\\AMRD","\\E-MMRD","\\G-MMRD",
//"","","\\EV","\\GM","\\E-MSW","\\G-MSW","\\E-MSD","\\G-MSD","\\E-AMRW","\\G-AMRW","\\E-AMRD","\\G-AMRD","\\MMRWW",
//"\\AMRWW","\\E-MMRWW","\\G-MMRWW","\\E-AMRWW","\\G-AMRWW","\\E-MSWW","\\G-MSWW"}; //各手法のファイルの名前
//// 使用する推定法のIDリスト
//int method_id_list[Method_N] = { 0,1,2,3,4,5,6,7,8,9,10,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30 };


//char file_method[16][64] = { "E-AMRw", "E-MMRw", "G-AMRw", "G-MMRw", "eAMRw", "eAMRwc", "eMMRw", "eMMRwc", "gAMRw", "gAMRwc", "gMMRw", "gMMRwc", "lAMRw", "lAMRwc", "lMMRw", "lMMRwc" };

void BubSort_Ascend(double x[][2], int item);										// ソート 昇順	Morm = 0;maximinのとき
void BubSort_Descend(double x[][2], int item);										//ソート　降順  Morm = 1;Maximaxのとき
void maximin(double u[Alt][N], double wR[N], double wL[N], double totalU[Alt], double z[Alt][N], int perm[Alt][N], int star[Alt]);	// 最小総合効用値の計算
int combi(int n, int m);														// ペアの番号
void arc_combi(int comb, int n[2]);

int main(void) {
	char max_or_min[2][64] = {"maximin","Maximax"};
	//maximin,Maximaxを定義する
	for (int Morm = 0;Morm < 2;Morm++) {
		//u1 or u2の二択
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


			FILE* fp_uM, * fp_iw;	// 入力ファイル

			char filepath[256];

			// 真の区間重要度setting
			for (rt = 0; rt < Rt; rt++) {

				printf("u%d:%s\n", utility, true_weight_list[rt]);

				// 真の区間重要度読み込み
				sprintf(filepath, "..\\data\\true_interval_weight_set%s\\%s\\Given_interval_weight.csv", eval_num_filepath[evaluation_num], true_weight_list[rt]);
				fp_iw = fopen(filepath, "r");
				if (fp_iw == nullptr) {
					printf("エラー: ファイルを開けませんでした (fp_iw): %s\n", filepath);
				}
				printf("開けた");
				for (f = 0; f < 2 * N; f++) {
					fscanf(fp_iw, "%lf,", &trueW[f]);
				}
				for (f = 0; f < 2 * N; f++) {
					if (f % 2 == 0) true_wL[f / 2] = trueW[f];
					else true_wR[f / 2] = trueW[f];
				}
				fclose(fp_iw);

				// 周辺効用値u読み込み
				sprintf(filepath, "..\\data\\効用値行列%s%s\\u.csv", which_utility[utility_num], util_num_filepath[evaluation_num]);
				fp_uM = fopen(filepath, "r");
				if (fp_uM == nullptr) {
					printf("エラー: ファイルを開けませんでした (fp_iw): %s\n", filepath);
				}
				printf("開けた");
				for (f = 0; f < Alt * M; f++) {
					fscanf(fp_uM, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &ukari[f][0], &ukari[f][1], &ukari[f][2], &ukari[f][3], &ukari[f][4], &ukari[f][5], &ukari[f][6], &ukari[f][7]);
				}
				fclose(fp_uM);

#pragma omp parallel for            //並行処理するための文章そのためにこの中でつかう変数はこの中で定義する
				for (md = 0; md < Method_N; md++) {	// 各推定法のループ
					int thread_id = omp_get_thread_num();
					printf("Thread %d is processing iteration %d\n", thread_id, md);
					int i, j, l, n, m, rp;			//ループ用変数
					FILE* fp_cs, * fp_rc, * fp_tc;	// 出力ファイル
					FILE* fp_ew;
					int k[2];
					int method_id;					// 推定法のループカウンタ
					double u[Alt][N];				// 周辺効用値u[Alt][N]

					double wR[N], wL[N];				// 推定区間重要度
					double wkari[R][2 * N];			// 推定区間重要度
					double wsum;
					double yR[N], yL[N];			// Y = [t*wL, t*wR]
					double z[Alt][N];				// 総合効用値で用いる重要度　z in Y
					double yR2[N], yL2[N];			// Y' = [r*yL, r*yR]

					double totalU[Alt];				// 総合効用値
					double totalU2[Alt];			// 総合効用値'

					double t_snap;					// t折れ曲がる点
					double t[101][3];				// t保管用　交差も含む
					double tL, tR;
					double sL, sR;					// 真の区間重要度の場合
					double s_snap;					// s折れ曲がる点
					double sh[101][3];				// s保管用　交差も含む
					double alphaR, alphaL;

					int perm[Alt][N];				// 周辺効用値が何番目に大きい要素か
					double udp[N][2];				// ある代替案の各基準の効用と元々どの基準かを記録
					int star[Alt];					// yi^L < zi < yi^R となるiの値
					int star2[Alt];
					double Sl, S[100];
					double r, r2[100][2];
					int it_t, it_s;
					int count1;
					int flag1;
					int flag2;						// 誤差判定
					double flag3[Alt];
					double swap1;
					int swap2;
					int rank1[100][Alt], rank2[100][Alt];
					double rank1c_s[100], rank2c_t[100];	// 順序変化s,t
					int count_rank1, count_rank2;
					double tbud[Alt][Alt], bud[Alt][Alt];
					int countC10[100][100];

					double store1, store2;
					char method_name[64];

					int id;
					char bff[1024];
					char* pbff, * ebff;
					char filepath_o[256];
					// 推定法
					method_id =md;

					//// 推定区間重要度の読みみ
					sprintf(filepath_o, "..\\data\\Simp\\Simp%s\\a%d\\%s%s\\Simp.csv", eval_num_filepath[evaluation_num], Generated_method, true_weight_list[rt], method_name_list[method_id]);
					fp_ew = fopen(filepath_o, "r");
					if (fp_ew == nullptr) {
						printf("エラー: ファイルを開けませんでした (fp_iw): %s\n", filepath_o);
					}
					fgets(bff, 1024, fp_ew);
					fgets(bff, 1024, fp_ew);
					fgets(bff, 1024, fp_ew);
					// ここまでSimp.csvのヘッダー処理, 以降, 推定区間重要度を読み込み
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


					// 出力ファイル
					sprintf(filepath_o, "..\\data\\a%d\\%s%s%s\\%s%s\\%s_count_pairs_in_squares_%d_u%d.csv", Generated_method, max_or_min[Morm], which_utility[utility_num], eval_num_filepath[evaluation_num], true_weight_list[rt], method_name_list[method_id], max_or_min[Morm],R,utility);
					fp_cs = fopen(filepath_o, "w");

					/*
					sprintf(filepath_o, "%s%s", method_name, "\\ranking_change.csv");
					fp_rc = fopen(filepath_o, "w");

					sprintf(filepath_o, "%s%s", method_name, "\\t_change.csv");
					fp_tc = fopen(filepath_o, "w");
					*/
					for (m = 0; m < M; m++) {

						// 周辺効用値の呼び出し
						for (i = 0; i < Alt; i++) {
							for (j = 0; j < N; j++) {
								u[i][j] = ukari[Alt * m + i][j];
							}
						}

						// 周辺効用値昇順
						for (j = 0; j < Alt; j++) {
							for (i = 0; i < N; i++) {
								udp[i][0] = u[j][i];	// iはどの基準かを表す
								udp[i][1] = i;
							}
							Morm == 0 ? BubSort_Ascend(udp, N) : BubSort_Descend(udp, N);
							for (i = 0; i < N; i++) {
								perm[j][i] = (int)udp[i][1];
							}
						}

						////真の区間重要度において正規性を満たす範囲sL,sRの計算
						alphaR = 10, alphaL = 0;	// 初期化
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
						// 初期化
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

							// マキシミン総合効用値
							maximin(u, yR, yL, totalU, z, perm, star);

							if (it_s == 1) { // 1回目だけrank1を計算
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

							// 次に傾きが変わる点を調べる
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

							// 微小誤差スキップ
							if (r == 1) {
								r = r - epsi;
								flag2 = -1;
							}

							// tLより小さくならないように
							if (r * s_snap < sL) {
								r = sL / s_snap;
								// printf("end\n");
							}

							// 総合効用値U'を計算
							for (i = 0; i < N; i++) {
								yR2[i] = yR[i] * r;
								yL2[i] = yL[i] * r;
							}

							maximin(u, yR2, yL2, totalU2, z, perm, star2);
							//printf("U2	:%f, %f, %f, %f, %f\n", totalU2[0], totalU2[1], totalU2[2], totalU2[3], totalU2[4]);

							// 入れ替わる点の有無 全ペア調査
							// 初期化
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

							// r'昇順
							if (count1 > 0) {
								BubSort_Ascend(r2, count1);
							}

							for (i = count1 - 1; i >= 0; i--) {
								sh[it_s][0] = s_snap * r2[i][0];
								sh[it_s][1] = 1;		// 交差
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
								sh[it_s][1] = flag2;			// 折れ
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
								sh[it_s][1] = -1;			// 誤差スキップ
								sh[it_s][2] = flag1;
								it_s++;
								printf("s:s=sr\n");
								flag2 = 1;	// 誤差スキップ
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
						// 出力
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

						// 推定区間重要度のループ
						for (rp = 0; rp < R; rp++) {

							// 推定区間重要度の呼び出し
							for (i = 0; i < N; i++) {
								wL[i] = wkari[rp][2 * i];
								wR[i] = wkari[rp][2 * i + 1];
							}

							// 推定区間重要度において正規性を満たす範囲tL,tRの計算
							alphaR = 10, alphaL = 0;	// 初期化
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
							// 初期化
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

								// マキシミン総合効用値
								maximin(u, yR, yL, totalU, z, perm, star);

								if (it_t == 1) { // 1回目だけrank2を計算
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

								// 次に傾きが変わる点を調べる
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

								// 微小誤差スキップ
								if (r == 1) {
									r = r - epsi;
									flag2 = -1;
								}

								// tLより小さくならないように
								if (r * t_snap < tL) {
									r = tL / t_snap;
									// printf("end\n");
								}

								// 総合効用値U'を計算
								for (i = 0; i < N; i++) {
									yR2[i] = yR[i] * r;
									yL2[i] = yL[i] * r;
								}

								maximin(u, yR2, yL2, totalU2, z, perm, star2);
								//printf("U2	:%f, %f, %f, %f, %f\n", totalU2[0], totalU2[1], totalU2[2], totalU2[3], totalU2[4]);

								// 入れ替わる点の有無 全ペア調査
								// 初期化
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

								// r'昇順
								if (count1 > 0) {
									BubSort_Ascend(r2, count1);
								}

								for (i = count1 - 1; i >= 0; i--) {
									t[it_t][0] = t_snap * r2[i][0];
									t[it_t][1] = 1;		// 交差
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
									t[it_t][1] = flag2;			// 折れ
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
									t[it_t][1] = flag2;			// 誤差スキップ
									t[it_t][2] = flag1;
									it_t++;
									//printf("t:t=tr\n");
									flag2 = 1;	// 誤差スキップ
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
							// 出力
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

							// 大小関係のカウント
							for (i = count_rank1 - 1; i >= 0; i--) {
								//真
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

								// 推定
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

/* maximin:前の要素の方が大きければ交換する*/
void BubSort_Ascend(double x[][2], int item) // 昇順ソート
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
/* maximax:前の要素の方がちいさければ交換する*/
void BubSort_Descend(double x[][2], int item) // 降順ソート
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
	double cap;				// 変数である重要度の総和が1となるようにするため

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

		//最適値の計算
		totalU[k] = 0;
		for (i = 0;i < N;i++) {
			totalU[k] = totalU[k] + u[k][i] * z[k][i];
		}
	}
}

int combi(int n, int m)
{
	/*
	例えば：代替案4個でペアの番号を振るとき
	0&1=0, 0&2=1, 0&3=2, 1&2=3, 1&3=4, 2&3=5　というふうに振れるように
	(Alt) C 2通り
	*/

	int i, o, p;
	int comb = 0;		// 組み合わせの番号

	if (n < m) {
		o = n + 1; p = m + 1;
	}
	else if (n > m) {
		o = m + 1; p = n + 1;
	}
	else {
		return -1;	// n,mが等しい場合はエラー
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