function summary_SDsynccells_frame = func_countSDsync_221105_2(X)
% X:２値化したデータで各frameに対する同期細胞数の行列

binarySD_SyncCount_Total = X>=2;%2cellの同期したframeを1, 1cell以下を0に
Sum_binarySD_SyncCount_Total = sum(binarySD_SyncCount_Total);%2cell以上の同期したframe数の合計

binarySD_SyncCount2 = X==2;%2cellの同期したframeを1, 1cell以下を0に
Sum_binarySD_SyncCount2 = sum(binarySD_SyncCount2);%2cell以上の同期したframe数の合計
binarySD_SyncCount3 = X==3;%3cellの同期したframeを1, 1cell以下を0に
Sum_binarySD_SyncCount3 = sum(binarySD_SyncCount3);%3cell以上の同期したframe数の合計
binarySD_SyncCount4 = X==4;%4cellの同期したframeを1, 1cell以下を0に
Sum_binarySD_SyncCount4 = sum(binarySD_SyncCount4);%4cell以上の同期したframe数の合計
binarySD_SyncCount5 = X>=5;%5cellの同期したframeを1, 1cell以下を0に
Sum_binarySD_SyncCount5 = sum(binarySD_SyncCount5);%5cell以上の同期したframe数の合計

binarySD_SyncCount1 = X==1;%2cellの同期したframeを1, 1cell以下を0に
Sum_binarySD_SyncCount1 = sum(binarySD_SyncCount1);%2cell以上の同期したframe数の合計
binarySD_Event = X>=1;%2cellの同期したframeを1, 1cell以下を0に
Sum_binarySD_Event = sum(binarySD_Event);%2cell以上の同期したframe数の合計


Num_SDsyncCells = ["2 cells" "3 cells" "4 cells" ">= 5 cells" "Total Sync" "1 cell" "Total event"];
Summary_binarySD_SyncCount = [Sum_binarySD_SyncCount2',Sum_binarySD_SyncCount3',Sum_binarySD_SyncCount4',Sum_binarySD_SyncCount5',Sum_binarySD_SyncCount_Total',Sum_binarySD_SyncCount1',Sum_binarySD_Event'];
Res_binarySD_SyncCount = [Num_SDsyncCells;Summary_binarySD_SyncCount];
summary_SDsynccells_frame = Res_binarySD_SyncCount;
