function nice_groupboxplot(data,lb,ub,names,group_name, gap, subtitle,label_vertical,theme_number)
%nice_groupboxplot( [rand(100,1) rand(100,1)+.2 rand(100,1)+.3 rand(100,1) rand(100,1)+.2 rand(100,1)+.3], -2,2,{'a','b','c'}, {'snr=1','snr=2'}, 2, '', '',1)
Y = data;
[~, entry_number] = size(Y);
item_number = length(names);
group_number = entry_number/item_number;

%color list
C1=[59 125 183;244 146 121;242 166 31;180 68 108;220 211 30]./255;
C2=[102,173,194;36,59,66;232,69,69;194,148,102;54,43,33]./255;
C3=[38,140,209;219,51,46;41,161,153;181,138,0;107,112,196]./255;
C4=[110,153,89;230,201,41;79,79,54;245,245,245;199,204,158]./255;
C5=[235,75,55;77,186,216;2,162,136;58,84,141;245,155,122]./255;
C6=[23,23,23;121,17,36;44,9,75;31,80,91;61,36,42]./255;
C7=[126,15,4;122,117,119;255,163,25;135,146,73;30,93,134]./255;

colorList=C4;

% 绘图
% y = Y(:);
% groups = ones(sample_number,1)*(1:entry_number);
% groups = groups(:);
% gap = 0.5/item_number;

%
% positions = [];
for i = 1:item_number
    %     positions = [positions i:gap:i+gap*(item_number-1)];
    datacell{i} = Y(:,i:item_number:(group_number-1)*item_number+i);
    names_empty{i} = '';
end



boxplotGroup(datacell, 'PrimaryLabels', names_empty, ...
    'SecondaryLabels',group_name, 'InterGroupSpace', gap, 'Symbol','o','OutlierSize',3,'Colors',[0,0,0]);
%boxplot(y,groups,'Symbol','o','OutlierSize',3,'Colors',[0,0,0], 'positions', positions);

% 坐标区域属性设置
ax=gca;hold on;
ax.LineWidth=1.1;
ax.FontSize=12;
ax.FontName='Arial';
ax.XTickLabel=names_empty;
ax.Title.String=subtitle;
ax.Title.FontSize=13;
ax.YLabel.String=label_vertical;

% 修改线条粗细
lineObj=findobj(gca,'Type','Line');
for i=1:length(lineObj)
    lineObj(i).LineWidth=1;
    lineObj(i).MarkerFaceColor=[1,1,1].*.3;
    lineObj(i).MarkerEdgeColor=[1,1,1].*.3;
end

% 为箱线图的框上色
boxObj=findobj(gca,'Tag','Box');
for i=1:length(boxObj)
    patch(boxObj(i).XData,boxObj(i).YData,colorList(ceil(i/group_number),:),'FaceAlpha',0.5,...
        'LineWidth',1.1);
end

%outliers labeling
% D = [];
% for i = 1:item_number
%     %     positions = [positions i:gap:i+gap*(item_number-1)];
%     D = [D Y(:,i:item_number:(group_number-1)*item_number+i)];
% end

columns = size(Y,2);
axis([0 1.5*entry_number-1.5 lb ub]);
ax = gca;
ax.YGrid = 'on';
ylim1=get(gca,'ylim');
% set(gca,'YTick',93:0.5:100)
ylabel(label_vertical)
n_outliers = sum(Y < lb);
for i = 1 : columns
    temp_outliers = n_outliers(i);
    if temp_outliers ~= 0
        str_temp = ['+', int2str(temp_outliers)];
        yvalue = lb + diff(ylim1)*.05;
        text(i+floor((i-1)/item_number)*gap + 0.1,yvalue,str_temp);
    end
end

% set theme
selfGrootDefault(theme_number)

names_flip = flip(names);
for i = 1:entry_number
    if mod(i,group_number) == 0
        lg{i} = names_flip{ceil(i/group_number)};
    else
        lg{i} = '';
    end
end

legend(lg)
end

function selfGrootDefault(theme)
if nargin<1
    theme=0;
end

try
    set(groot,'defaultLineLineWidth',2)
    set(groot,'defaultAxesFontName','Arial')
    set(groot,'defaultAxesColorOrder',[102,194,165; ...
        252,140,98; ...
        142,160,204; ...
        231,138,195; ...
        166,217,83; ...
        255,217,48; ...
        229,196,148; ...
        179,179,179]./255)
    % set(groot,'defaultAxesProjection','perspective')
    set(groot,'defaultSurfaceEdgeColor',[1,1,1].*.7);
    set(groot,'defaultSurfaceEdgeAlpha',.5);


    CM=[127 0 255 126 0 254 125 1 254 125 2 254 124 3 254 123 5 254 122 6 254 121 7 254 121 8 254 120 10 254 119 11 254 118 12 254 117 13 254 117 14 254 116 16 254 115 17 254 114 18 254 113 19 254 113 20 254 112 22 254 111 24 254 110 25 254 109 26 254 108 28 254 108 29 254 107 30 254 106 31 254 105 32 254 104 34 254 104 35 254 103 36 254 102 37 254 101 38 254 100 40 254 100 41 254 99 42 254 98 43 253 97 45 253 96 46 253 96 47 253 95 48 253 94 50 253 93 51 253 92 53 253 92 54 253 91 55 253 90 56 253 89 58 253 88 59 253 87 60 253 87 61 253 86 62 252 85 64 252 84 65 252 83 66 252 83 67 252 82 69 252 81 70 252 80 71 252 79 72 252 79 73 252 78 75 251 77 76 251 76 77 251 75 78 251 75 79 251 74 81 251 73 82 251 72 83 251 71 84 251 70 86 250 70 87 250 69 88 250 68 89 250 67 90 250 66 92 250 66 93 250 65 94 250 64 95 250
        63 96 250 62 97 249 62 98 249 61 99 249 60 100 249 59 102 249 58 103 249 58 104 249 57 105 249 56 106 248 55 108 248 54 109 248 53 110 248 53 111 248 52 112 248 51 114 248 50 115 247 49 116 247 49 116 247 48 118 247 47 119 247 46 120 246 45 121 246 45 122 246 44 124 246 43 125 246 42 126 246 41 127 246 41 127 246 40 129 245 39 130 245 38 131 245 37 132 245 36 134 244 36 134 244 35 135 244 34 136 244 33 137 244 32 139 243 32 140 243 31 141 243 30 142 243 29 143 243 28 144 243 28 145 243 27 146 243 26 147 242 25 148 242 24 149 242 24 150 242 23 151 242 22 152 241 21 153 241 20 154 241 20 155 241 19 156 241 18 157 240 17 158 240 16 159 239 15 160 239 15 161 239 14 162 239 13 163 239 12 164 238 11 165 238 11 166 238 10 167 238 9 168 238 8 169 237 7 170 237 7 170 237 6 172 237 5 173 237 4 174 236 3 175 236 3 175 236 2 176 236 1 177 236 0 178 235
        0 179 235 0 180 234 0 181 234 1 182 234 2 183 234 3 184 234 4 185 233 4 185 233 5 186 233 6 187 232 7 188 232 8 189 232 8 189 232 9 190 232 10 191 231 11 192 231 12 193 230 12 193 230 13 194 230 14 195 230 15 196 230 16 197 229 17 198 229 17 198 229 18 199 228 19 200 228 20 201 228 21 202 228 21 202 228 22 203 227 23 204 227 24 205 226 25 206 226 25 206 226 26 207 226 27 208 226 28 209 225 29 209 225 29 209 225 30 210 224 31 211 224 32 212 223 33 213 223 33 213 223 34 214 223 35 214 223 36 215 222 37 216 222 38 217 221 38 217 221 39 218 221 40 219 220 41 219 220 42 220 220 42 220 220 43 221 220 44 222 219 45 222 219 46 223 218 46 223 218 47 224 218 48 225 217 49 225 217 50 226 216 50 226 216 51 227 216 52 228 215 53 228 215 54 229 215 55 229 215 55 229 215 56 230 214 57 231 214 58 232 213 59 232 213 59 232 213 60 233 212 61 233 212 62 234 211
        63 235 211 63 235 211 64 236 210 65 236 210 66 237 209 67 237 209 67 237 209 68 238 209 69 238 209 70 239 208 71 239 208 72 240 207 72 240 207 73 240 207 74 241 206 75 241 206 76 242 205 76 242 205 77 242 205 78 243 204 79 243 204 80 244 203 80 244 203 81 244 203 82 245 202 83 245 202 84 246 201 84 246 201 85 246 201 86 246 200 87 246 200 88 247 199 88 247 199 89 247 199 90 248 198 91 248 198 92 249 197 93 249 197 93 249 197 94 249 196 95 249 196 96 250 195 97 250 195 97 250 195 98 250 194 99 250 194 100 251 193 101 251 193 101 251 193 102 251 192 103 251 192 104 252 191 105 252 191 105 252 191 106 252 190 107 252 190 108 253 189 109 253 189 110 253 188 110 253 188 111 253 188 112 253 187 113 253 187 114 254 186 114 254 186 115 254 186 116 254 185 117 254 185 118 254 184 118 254 184 119 254 184 120 254 183 121 254 183 122 254 182 122 254 182 123 254 181 124 254 180 125 254 180 126 254 179
        127 254 179 127 254 179 128 254 178 129 254 178 130 254 177 131 254 177 131 254 177 132 254 176 133 254 176 134 254 175 135 254 175 135 254 175 136 254 174 137 254 174 138 254 173 139 254 172 139 254 172 140 253 171 141 253 171 142 253 170 143 253 170 143 253 170 144 253 169 145 253 169 146 252 168 147 252 168 148 252 167 148 252 167 149 252 167 150 251 166 151 251 165 152 251 164 152 251 164 153 251 164 154 250 163 155 250 163 156 250 162 156 250 162 157 250 162 158 249 161 159 249 161 160 249 160 160 249 160 161 249 159 162 248 158 163 248 158 164 247 157 165 247 157 165 247 157 166 246 156 167 246 156 168 246 155 169 246 154 169 246 154 170 245 153 171 245 153 172 244 152 173 244 152 173 244 152 174 243 151 175 243 151 176 242 150 177 242 149 177 242 149 178 241 148 179 241 148 180 240 147 181 240 147 182 239 146 182 239 146 183 239 146 184 238 145 185 238 144 186 237 143 186 237 143 187 237 143 188 236 142 189 236 142 190 235 141
        190 235 141 191 234 140 192 233 139 193 233 139 194 232 138 194 232 138 195 232 138 196 231 137 197 230 136 198 229 135 199 229 135 199 229 135 200 228 134 201 228 134 202 227 133 203 226 132 203 226 132 204 225 131 205 225 131 206 224 130 207 223 130 207 223 130 208 222 129 209 222 128 210 221 127 211 220 127 211 220 127 212 219 126 213 219 126 214 218 125 215 217 124 215 217 124 216 216 123 217 215 123 218 214 122 219 214 122 220 213 121 220 213 121 221 212 120 222 211 119 223 210 119 224 209 118 224 209 118 225 209 117 226 208 116 227 207 116 228 206 115 228 206 115 229 205 115 230 204 114 231 203 113 232 202 112 232 202 112 233 201 112 234 200 111 235 199 110 236 198 109 237 197 109 237 197 109 238 196 108 239 195 108 240 194 107 241 193 106 241 193 106 242 192 105 243 191 105 244 190 104 245 189 103 245 189 103 246 188 102 247 187 102 248 186 101 249 185 100 249 185 100 250 184 99 251 183 99 252 182 98 253 181 97 253 180 97
        254 179 96 254 178 96 255 177 95 255 176 95 255 175 94 255 175 94 255 174 93 255 173 92 255 172 92 255 170 91 255 170 91 255 169 90 255 168 89 255 167 89 255 166 88 255 165 88 255 164 87 255 163 86 255 162 86 255 161 85 255 160 85 255 159 84 255 158 83 255 157 83 255 156 82 255 155 81 255 154 81 255 153 80 255 152 80 255 151 79 255 150 78 255 149 78 255 148 77 255 147 77 255 146 76 255 145 75 255 144 75 255 143 74 255 142 74 255 141 73 255 140 72 255 139 72 255 137 71 255 136 71 255 135 70 255 134 69 255 133 68 255 132 68 255 131 68 255 130 67 255 129 66 255 127 65 255 127 65 255 126 65 255 125 64 255 124 63 255 122 62 255 121 62 255 120 62 255 119 61 255 117 60 255 116 59 255 116 59 255 115 59 255 114 58 255 112 57 255 111 56 255 110 56 255 109 56 255 108 55 255 106 54 255 105 53 255 104 53 255 103 53 255 102 52 255 100 51 255 99 50 255 98 49 255 97 49
        255 96 48 255 95 48 255 94 47 255 93 46 255 92 46 255 90 45 255 89 45 255 88 44 255 87 43 255 86 43 255 84 42 255 83 42 255 82 41 255 81 40 255 79 39 255 78 39 255 77 39 255 76 38 255 75 37 255 73 36 255 72 36 255 71 36 255 70 35 255 69 34 255 67 33 255 66 33 255 65 32 255 64 31 255 62 31 255 61 30 255 60 30 255 59 29 255 58 28 255 56 28 255 55 27 255 54 26 255 53 26 255 51 25 255 50 25 255 48 24 255 47 23 255 46 23 255 44 22 255 43 21 255 42 20 255 41 20 255 40 20 255 38 19 255 37 18 255 36 17 255 35 17 255 34 17 255 32 16 255 31 15 255 30 14 255 29 14 255 28 14 255 26 13 255 25 12 255 24 11 255 22 10 255 20 9 255 19 9 255 18 9 255 17 8 255 16 7 255 14 6 255 13 6 255 12 6 255 11 5 255 10 4 255 8 3 255 7 3 255 6 3 255 5 2 255 3 1 255 2 0 255 1 0 255 0 0]./255;
    set(groot,'defaultFigureColormap',reshape(CM',3,[])');
    set(groot,'defaultPatchLineWidth',2);
    set(groot,'defaultPatchFaceAlpha',.5);
    set(groot,'defaultAreaFaceAlpha',.6)
    set(groot,'defaultAreaLineWidth',1.5)
    set(groot,'defaultAreaEdgeColor',[.2,.2,.2])
    set(groot,'defaultBarLineWidth',1.5)
    set(groot,'defaultBarFaceAlpha',.8)
    set(groot,'defaultStemLineWidth',1.2)
catch
end


oriDefault()
switch theme
    case {0,'origin'}
        set(groot,'defaultLineLineWidth',0.5)
        set(groot,'defaultAxesColorOrder',[    0    0.4470    0.7410;...
            0.8500    0.3250    0.0980;...
            0.9290    0.6940    0.1250;...
            0.4940    0.1840    0.5560;...
            0.4660    0.6740    0.1880;...
            0.3010    0.7450    0.9330;...
            0.6350    0.0780    0.1840])
        set(groot,'defaultLegendTextColor',[0,0,0])
        set(groot,'defaultAxesProjection','orthographic')
        set(groot,'defaultFigureColormap',parula)
        set(groot,'defaultSurfaceEdgeColor',[0,0,0]);
        set(groot,'defaultSurfaceEdgeAlpha',1);
        set(groot,'defaultPatchLineWidth',.5);
        set(groot,'defaultPatchFaceAlpha',1);
        set(groot,'defaultAxesFontSize',10)
        set(groot,'defaultAreaFaceAlpha',1)
        set(groot,'defaultAreaLineWidth',.5)
        set(groot,'defaultAreaEdgeColor',[0,0,0])
        set(groot,'defaultBarLineWidth',.5)
        set(groot,'defaultBarFaceAlpha',1)
        set(groot,'defaultStemLineWidth',.5)
    case {1,'gbase'}
        set(groot,'defaultAxesLineWidth',1.2)
        set(groot,'defaultAxesXMinorTick','on')
        set(groot,'defaultAxesYMinorTick','on')
        set(groot,'defaultAxesZMinorTick','on')
        set(groot,'defaultAxesGridLineStyle','-.')
        set(groot,'defaultAxesXColor',[1,1,1].*.3)
        set(groot,'defaultAxesYColor',[1,1,1].*.3)
        set(groot,'defaultAxesZColor',[1,1,1].*.3)
        set(groot,'defaultAxesXGrid','on')
        set(groot,'defaultAxesYGrid','on')
        set(groot,'defaultAxesZGrid','on')
        set(groot,'defaultAxesTickDir','in','defaultAxesTickDirMode','auto')
        set(groot,'defaultLegendTextColor',[0,0,0])
    case {2,'gbase2'}
        set(groot,'defaultAxesLineWidth',1.2)
        set(groot,'defaultAxesXMinorTick','on')
        set(groot,'defaultAxesYMinorTick','on')
        set(groot,'defaultAxesZMinorTick','on')
        set(groot,'defaultAxesGridLineStyle','-.')
        set(groot,'defaultAxesXColor',[1,1,1].*.3)
        set(groot,'defaultAxesYColor',[1,1,1].*.3)
        set(groot,'defaultAxesZColor',[1,1,1].*.3)
        set(groot,'defaultAxesXGrid','on')
        set(groot,'defaultAxesYGrid','on')
        set(groot,'defaultAxesZGrid','on')
        set(groot,'defaultAxesColor',[249,250,245]./255)
        set(groot,'defaultAxesTickDir','in','defaultAxesTickDirMode','auto')
        set(groot,'defaultLegendTextColor',[0,0,0])
    case {3,'dark'}
        set(groot,'defaultAxesColor',[0,0,0])
        set(groot,'defaultFigureColor',[0,0,0])
        set(groot,'defaultAxesXMinorTick','on')
        set(groot,'defaultAxesYMinorTick','on')
        set(groot,'defaultAxesZMinorTick','on')
        set(groot,'defaultAxesXGrid','on')
        set(groot,'defaultAxesYGrid','on')
        set(groot,'defaultAxesZGrid','on')
        set(groot,'defaultAxesXMinorGrid','on','defaultAxesXMinorGridMode','manual')
        set(groot,'defaultAxesYMinorGrid','on','defaultAxesYMinorGridMode','manual')
        set(groot,'defaultAxesZMinorGrid','on','defaultAxesZMinorGridMode','manual')
        set(groot,'defaultFigureInvertHardCopy','off')
        set(groot,'defaultAxesGridColor',[0.2353 0.4235 0.4745])
        set(groot,'defaultAxesMinorGridColor',[0.1569 0.2039 0.1882])
        set(groot,'defaultAxesXColor',[1,1,1])
        set(groot,'defaultAxesYColor',[1,1,1])
        set(groot,'defaultAxesZColor',[1,1,1])
        set(groot,'defaultAxesLineWidth',1.2)
        set(groot,'defaultAxesGridAlpha',.8)
        set(groot,'defaultAxesMinorGridAlpha',0.7)
        set(groot,'defaultAxesTickDir','in','defaultAxesTickDirMode','auto')
        set(groot,'defaultLegendTextColor',[1,1,1])
    case {4,'dark2'}
        set(groot,'defaultAxesColor',[2,34,57]./255)
        set(groot,'defaultFigureColor',[2,34,57]./255)
        set(groot,'defaultFigureInvertHardCopy','off')
        set(groot,'defaultAxesXColor',[144,164,178]./255)
        set(groot,'defaultAxesYColor',[144,164,178]./255)
        set(groot,'defaultAxesZColor',[144,164,178]./255)
        set(groot,'defaultAxesLineWidth',1.2)
        set(groot,'defaultAxesXGrid','on')
        set(groot,'defaultAxesYGrid','on')
        set(groot,'defaultAxesZGrid','on')
        set(groot,'defaultAxesGridAlpha',.2)
        set(groot,'defaultAxesGridColor',[144,164,178]./255)
        set(groot,'defaultAxesXMinorTick','on')
        set(groot,'defaultAxesYMinorTick','on')
        set(groot,'defaultAxesZMinorTick','on')
        set(groot,'defaultAxesGridLineStyle','-.')
        set(groot,'defaultAxesTickDir','in','defaultAxesTickDirMode','auto')
        set(groot,'defaultLegendTextColor',[1,1,1])
    case {5,'ggray2'}
        set(groot,'defaultAxesColor',[234,234,242]./255)
        set(groot,'defaultFigureColor',[1,1,1])
        set(groot,'defaultAxesTickDir','out','defaultAxesTickDirMode','manual')
        set(groot,'defaultAxesXGrid','on')
        set(groot,'defaultAxesYGrid','on')
        set(groot,'defaultAxesZGrid','on')
        set(groot,'defaultAxesGridColor',[1,1,1])
        set(groot,'defaultAxesGridAlpha',1)
        set(groot,'defaultAxesLineWidth',1.2)
        set(groot,'defaultAxesXColor',[1,1,1].*.33)
        set(groot,'defaultAxesYColor',[1,1,1].*.33)
        set(groot,'defaultAxesZColor',[1,1,1].*.33)
        set(groot,'defaultAxesTickLength',[.005,.025])
        set(groot,'defaultLegendTextColor',[0,0,0])
        set(groot,'defaultAxesBox','off')
    case {6,'economist'}
        set(groot,'defaultAxesColor',[0.8400 0.8900 0.9200])
        set(groot,'defaultFigureColor',[0.8400 0.8900 0.9200])
        set(groot,'defaultFigureInvertHardCopy','off')
        set(groot,'defaultAxesBox','off')
        set(groot,'defaultAxesXGrid','off')
        set(groot,'defaultAxesYGrid','on')
        set(groot,'defaultAxesZGrid','on')
        set(groot,'defaultAxesGridColor',[1,1,1])
        set(groot,'defaultAxesGridAlpha',1)
        set(groot,'defaultAxesLineWidth',1.2)
        set(groot,'defaultAxesXColor',[1,1,1].*.33)
        set(groot,'defaultAxesYColor',[1,1,1].*.33)
        set(groot,'defaultAxesZColor',[1,1,1].*.33)
        set(groot,'defaultAxesTickLength',[.005,.025])
    case {7,'wsj'}
        set(groot,'defaultAxesColor',[0.9700 0.9500 0.8900])
        set(groot,'defaultFigureColor',[0.9700 0.9500 0.8900])
        set(groot,'defaultFigureInvertHardCopy','off')
        set(groot,'defaultAxesBox','off')
        set(groot,'defaultAxesXGrid','off')
        set(groot,'defaultAxesYGrid','on')
        set(groot,'defaultAxesZGrid','on')
        set(groot,'defaultAxesLineWidth',1)
        set(groot,'defaultAxesGridAlpha',1)
        set(groot,'defaultAxesTickLength',[.005,.025])
        set(groot,'defaultAxesGridLineStyle',':')
end
    function oriDefault(~,~)
        set(groot,'defaultAxesLineWidth',0.5)
        set(groot,'defaultAxesXMinorTick','off')
        set(groot,'defaultAxesYMinorTick','off')
        set(groot,'defaultAxesZMinorTick','off')
        set(groot,'defaultAxesXMinorGrid','off','defaultAxesXMinorGridMode','auto')
        set(groot,'defaultAxesYMinorGrid','off','defaultAxesYMinorGridMode','auto')
        set(groot,'defaultAxesZMinorGrid','off','defaultAxesZMinorGridMode','auto')
        set(groot,'defaultAxesGridLineStyle','-')
        set(groot,'defaultAxesGridColor',[1,1,1].*.15)
        set(groot,'defaultAxesMinorGridColor',[1,1,1].*.1)
        set(groot,'defaultAxesXColor',[1,1,1].*.15)
        set(groot,'defaultAxesYColor',[1,1,1].*.15)
        set(groot,'defaultAxesZColor',[1,1,1].*.15)
        set(groot,'defaultAxesGridAlpha',.15)
        set(groot,'defaultAxesMinorGridAlpha',.25)
        set(groot,'defaultAxesFontName','Helvetica')
        set(groot,'defaultAxesXGrid','off')
        set(groot,'defaultAxesYGrid','off')
        set(groot,'defaultAxesColor',[1,1,1])
        set(groot,'defaultFigureColor',[1,1,1].*.94)
        set(groot,'defaultAxesBox','off')
        set(groot,'defaultFigureInvertHardCopy','on')
        set(groot,'defaultAxesTickDir','in','defaultAxesTickDirMode','auto')
        set(groot,'defaultAxesTickLength',[.01,.025])

    end
end

function handles = boxplotGroup(varargin)
% BOXPLOTGROUP groups boxplots together with horizontal space between groups.
%   boxplotGroup(x) receives a 1xm cell array where each element is a matrix with
%   n columns and produces n groups of boxplot boxes with m boxes per group.
%
%   boxplotGroup(ax,x,___) specifies the axis handle, otherwise current axis is used.
%
%   boxplotGroup(___,'interGroupSpace',d) separates groups by d units along the x axis
%   where d is a positive, scalar integer (default = 1)
%
%   boxplotGroup(___,'groupLines', true) adds vertical divider lines between groups
%   (requires >=r2018b).
%
%   boxplotGroup(___,'primaryLabels', c) specifies the x tick label for each boxplot.
%   c is a string array or cell array of characters and must have one element per
%   box or one element per group-member. When undefined or when c is an empty cell {},
%   the x-axis is labeled with default x-tick labels.
%
%   boxplotGroup(___,'secondaryLabels', s) specifies the group labels for the boxplot
%   groups.  s is a string array or cell array of characters and must have one element
%   per group (see 'groupLabelType'). Ignored when s is an empty cell {}.
%
%   boxplotGroup(___,'groupLabelType', str) specifies how to label the groups by one of
%   the following options.
%    * 'horizontal': Group labels will be centered under the primary labels using a 2nd
%       invisible axis underlying the main axis (not supported in uifigures). To remove
%       the primary labels and only show secondary labels, set primary labels to empty
%       cell-strings (e.g. {'','',''}) or strings without characters (e.g. ["" "" ""]).
%    * 'vertical': Group labels will be vertical, between groups (requires Matlab >=2018b)
%    * 'both': Both methods will be used.
%
%   boxplotGroup(___, 'PARAM1', val1, 'PARAM2, val2, ...) sends optional name/value pairs
%   to the boxplot() function. Accepted parameters are BoxStyle, Colors, MedianStyle,
%   Notch, OutlierSize, PlotStyle, Symbol, Widths, DataLim, ExtremeMode, Jitter, and Whisker.
%   See boxplot documentation for details.
%
%   boxplotGroup(___, 'Colors', ___, 'GroupType', type) determines how to apply
%   colors to the groups.  'Colors' is a property of boxplots (see boxplot documentation).
%   When the colors value specifies multiple colors, the 'GroupType' determines how
%   the colors are distributed based on the following two options.
%    * 'betweenGroups' assigns color n to the n^th boxplot within each group (default).
%    * 'withinGroups' assigns color n to all boxplots within the n^th group.
%
%   h = boxplotGroup(___) outputs a structure of graphics handles.
%
% NOTE: If you're working with a grouping variable 'g', use the syntax boxplot(x,g) along
%   with the "Group Appearance" options described in Matlab's boxplot() documentation.
%   https://www.mathworks.com/help/stats/boxplot.html#d118e146984
%
% EXAMPLES:
% data = {rand(100,4), rand(20,4)*.8, rand(1000,4)*1.2};
%
% Required inputs
%   boxplotGroup(data)
%
% Set space between groups
%   boxplotGroup(data, 'interGroupSpace', 3)
%
% Specify labels and draw divider line
%   boxplotGroup(data, 'groupLines', true, 'PrimaryLabels', {'a' 'b' 'c'},...
%       'SecondaryLabels', {'Lancaster', 'Cincinnati', 'Sofia', 'Rochester'})
%
% Label groups with vertical lables
%   boxplotGroup(data, 'PrimaryLabels', {'a' 'b' 'c'}, 'SecondaryLabels', ...
%       {'Lancaster', 'Cincinnati', 'Sofia', 'Rochester'}, 'groupLabelType', 'vertical')
%
% Pass additional boxplot properties
%   boxplotGroup(data, 'PrimaryLabels', {'a' 'b' 'c'}, 'SecondaryLabels', ...
%       {'Lancaster', 'Cincinnati', 'Sofia', 'Rochester'}, 'groupLabelType', 'vertical', ...
%       'BoxStyle', 'filled', 'PlotStyle', 'Compact')
%
%
% Contact adam.danz@gmail.com for questions, bugs, suggestions, and high-fives.
% Copyright (c) 2020, Adam Danz  adam.danz@gmail.com
% All rights reserved
% Source: https://www.mathworks.com/matlabcentral/fileexchange/74437-boxplotgroup
% Changes history
% 200306 - v1.0.0 first release.
% 200308 - v1.1.0 Added recommendation to use boxplot() with grouping variable.
%                 Added axis handle as input to boxplot() call. Linkaxes changed
%                 from x to xy. Set axis2.Units to axis.Units.  Using linkprop
%                 to link position etc of main axis and axis2. Added DeleteFcn
%                 to main axis. Disabled toolbar for axis2. Added listener to
%                 resize axis2 when main axis is resized. Changes to help section.
% 200309 - v1.2.0 When 2nd axis is added, main axis is set to current axis.
% 200309 - v1.2.1 Suppress linkprops() and changes to toolbar suppression to work
%                 with versions prior to r2018b.
% 200309 - v1.2.2 Instead of creating new axis, default axis is gca().
% 210427 - v2.0.0 oncleanup returns hold state instead of conditional.  Added GroupType
%                 option and colorexpansion. Suppresses output unless requested.
%                 Checks matlab vs with xline(). Removed listener, storing hlink in axis.
%                 boxplot name-val arg check. Removing boxplot placeholders. XTicks now auto
%                 if labels aren't provided. Outputs now include boxplotGroup; vertical
%                 labels now the same fontsize and weight as axis font; Primary and secondary
%                 labels can be empty cell to ignore. Secondary labels now match ax1 font size,
%                 weight and name.
%% Check for axis handle in first input
if ~isempty(varargin) && ~isempty(varargin{1}) && isobject(varargin{1}(1)) % [3]
    if isgraphics(varargin{1}(1), 'axes')
        % first input is an axis
        h.axis = varargin{1} ;
        varargin(1) = [];
    else
        error('MATLAB:hg:InvalidHandle', 'Invalid handle')
    end
else
    h.axis = [];
end
%% Parse inputs
p = inputParser();
p.FunctionName = mfilename;
p.KeepUnmatched = true;	%accept additional parameter value inputs (passed to boxplot())
addRequired(p, 'x', @(x)validateattributes(x,{'cell'},{'row','nonempty'}))
addParameter(p, 'interGroupSpace', 1, @(x)validateattributes(x,{'double'},{'scalar','integer'}))
addParameter(p, 'primarylabels', [], @(x)validateattributes(x,{'string','cell'},{}))
addParameter(p, 'secondarylabels', [], @(x)validateattributes(x,{'string','cell'},{}))
addParameter(p, 'groupLines', false, @(x)validateattributes(x,{'logical','double'},{'binary'}))
addParameter(p, 'groupLabelType', 'Horizontal', @(x)ischar(validatestring(lower(x),{'vertical','horizontal','both'})))
addParameter(p, 'GroupType', 'betweenGroups', @(x)ischar(validatestring(lower(x),{'betweengroups','withingroups'})))
parse(p,varargin{:})
% Prepare the unmatched boxplot() parameters.
% If a param is passed that isn't accepted by boxplot(), an error is thrown from boxplot() function.
unmatchNameVal = reshape([fieldnames(p.Unmatched)'; struct2cell(p.Unmatched)'], 1, []);
% Check boxplot name-value parameters; group params, Position, and labels are not accepted.
supportedParams = {'BoxStyle','Colors','MedianStyle','Notch','OutlierSize','PlotStyle','Symbol','Widths', ...
    'DataLim','ExtremeMode','Jitter','Whisker'};
argOK = arrayfun(@(i)any(strncmpi(unmatchNameVal{i},supportedParams,numel(unmatchNameVal{i}))),...
    1:2:numel(unmatchNameVal)); % look for partial match
assert(all(argOK),'Parameter(s) not accepted in %s: [%s].', ...
    mfilename, strjoin(unmatchNameVal(find(~argOK)*2-1),', '))
% Check that each element of x is a matrix
assert(all(cellfun(@ismatrix, p.Results.x)), 'All elements of the cell array ''x'' must be a matrix.')
% Check that each matrix contains the same number of columns.
assert(numel(unique(cellfun(@(m)size(m,2),p.Results.x))) == 1, ...
    ['All elements of the cell array ''x'' must contain the same number of columns. '...
    'Pad the matricies that contain fewer columns with NaN values.']);
nargoutchk(0,1)
%% Compute horizontal spacing & check labels
nGroups = size(p.Results.x{1},2);       % number of columns of data / number of groups
nMembers = numel(p.Results.x);          % number of members per group
maxX = ((nMembers + p.Results.interGroupSpace) * nGroups) - p.Results.interGroupSpace;
xInterval = nMembers + p.Results.interGroupSpace;
% Check that labels (if any) are the right size
% PrimaryLabels: either 1 per group-member or 1 for each bar
if ~isempty(p.Results.primarylabels)
    assert(ismember(numel(p.Results.primarylabels),[nMembers, nMembers*nGroups]), ...
        sprintf(['The number of primary labels must equal either the number of bars per group (%d) '...
        'or the number of total bars (%d).'], nMembers, nMembers*nGroups))
end
% SecondaryLabels: 1 per group
if ~isempty(p.Results.secondarylabels)
    assert(isequal(numel(p.Results.secondarylabels),nGroups), ...
        sprintf('The number of secondary labels must equal either the number groups (%d).',nGroups))
end
% If all primary labels are empty chars do not add the newline to secondary labels.
if ~isempty(p.Results.primarylabels) &&  all(cellfun(@isempty,cellstr(p.Results.primarylabels)))
    horizSecondaryLabelAddon = '';
else
    horizSecondaryLabelAddon = '\newline';
end
%% Set colors
% Assumes ColorGroup property is not specified.
colorsIdx = strcmpi('Colors',unmatchNameVal);
if any(colorsIdx)
    cvalIdx = find(colorsIdx,1,'first')+1;
    if isempty(unmatchNameVal{cvalIdx})
        % Colors val is empty; remove Colors name-val pair
        unmatchNameVal(cvalIdx-[1,0]) = [];
    else
        unmatchNameVal{cvalIdx} = colorexpansion(unmatchNameVal{cvalIdx}, p, nGroups, nMembers);
    end
end
%% Do plotting
if isempty(h.axis)
    h.axis = gca();
end
h.figure = ancestor(h.axis,'figure');
isTiledLayout = strcmpi(h.axis.Parent.Type,'tiledlayout');
if isTiledLayout % [6]
    origTLOState = warning('query', 'MATLAB:handle_graphics:Layout:NoPositionSetInTiledChartLayout');
    TLOcleanup = onCleanup(@()warning(origTLOState));
    warning('off','MATLAB:handle_graphics:Layout:NoPositionSetInTiledChartLayout')
end
% Store pre-existing boxplot object handles
bptag = 'boxplot'; % tag Matlab assigns to bp group
bpobjPre = findobj(h.axis,'tag',bptag);
originalHoldStatus = ishold(h.axis);
holdStates = {'off','on'};
returnHoldState = onCleanup(@()hold(h.axis,holdStates{originalHoldStatus+1}));
hold(h.axis, 'on')
x = cell(1,nMembers);
existingTextObjs = findobj(h.axis,'Type','Text');
for i = 1:nMembers
    x{i} = i : xInterval : maxX;
    temp = nan(size(p.Results.x{i},1), max(x{i}));
    temp(:,x{i}) = p.Results.x{i};
    boxplot(h.axis, temp, unmatchNameVal{:})
end
% Remove dummy boxplots placeholders
bpobjNew = findobj(h.axis,'tag',bptag);
bpobjNew(ismember(bpobjNew, bpobjPre)) = [];
for g = 1:numel(bpobjNew)
    tags = unique(get(bpobjNew(g).Children,'tag'),'stable');
    tags(cellfun(@isempty,tags)) = [];
    for j = 1:numel(tags)
        obj = findobj(bpobjNew(g),'tag',tags{j});
        obj(~isprop(obj,'YData')) = [];
        YData = get(obj,'YData');
        if ~iscell(YData)
            YData = {YData};
        end
        isDummy = cellfun(@(c)all(isnan(c),2),YData);
        delete(obj(isDummy))
    end
end
axis(h.axis, 'tight')
limGap = (p.Results.interGroupSpace+1)/2;
set(h.axis,'XTickMode','Auto','XTickLabelMode','Auto','xlim',[1-limGap, maxX+limGap]) %[1]
yl = ylim(h.axis);
ylim(h.axis, yl + [-range(yl)*.05, range(yl)*.05])
% Remove boxplot's text-tics [1]
allTextObjs = findobj(h.axis,'Type','Text');
isBoxplotText = ~ismember(allTextObjs,existingTextObjs);
set(allTextObjs(isBoxplotText), 'String','','Visible','off')
% Set primary labels if provided
if ~isempty(p.Results.primarylabels)
    h.axis.XTick = sort([x{:}]);
    h.axis.XTickLabel = p.Results.primarylabels;
end
% Set secondary labels if provided
vertLinesDrawn = false;
groupLabelType = p.Results.groupLabelType;
if ~isempty(p.Results.secondarylabels)
    if any(strcmpi(groupLabelType, {'horizontal','both'}))
        % Try to detect figure type [4]
        if verLessThan('Matlab','9.0')      %version < 16a (release of uifigs)
            isuifig = @(~)false;
        elseif verLessThan('Matlab','9.5')  % 16a <= version < 18b
            isuifig = @(h)~isempty(matlab.ui.internal.dialog.DialogHelper.getFigureID(h));
        else                                % version >= 18b (written in r21a)
            isuifig = @(h)matlab.ui.internal.isUIFigure(h) && ~isprop(h,'LiveEditorRunTimeFigure');
        end
        isUIFigure = isuifig(h.figure);
        if isUIFigure
            groupLabelType = 'vertical';
            warning('BOXPLOTGRP:uifig','''Horizontal'' GroupLabelType is not supported with UIFIgures. GroupLabelType was changed to ''vertical''.')
        else
            % Tick label rotation must be 0 if using both primary & secondary horizontal labels
            h.axis.XAxis.TickLabelRotation = 0;
            % Compute x position of secondary labels
            if isa(h.axis,'matlab.ui.control.UIAxes')
                axFcn = @uiaxes;
            else
                axFcn = @axes;
            end
            if verLessThan('Matlab','9.8') %r2020a
                posProp = 'Position';
            else
                posProp = 'InnerPosition';
            end
            secondaryX = (nMembers : nMembers + p.Results.interGroupSpace : maxX) - (nMembers-1)/2;
            secondaryLabels = strcat(horizSecondaryLabelAddon,p.Results.secondarylabels); %[2]
            h.axis2 = axFcn(h.figure,'Units',h.axis.Units, 'OuterPosition', h.axis.OuterPosition, ...
                'ActivePositionProperty', h.axis.ActivePositionProperty,'xlim', h.axis.XLim, ...
                'TickLength', [0 0], 'ytick', [], 'Color', 'none', 'XTick', secondaryX, ...
                'TickLabelInterpreter','tex','XTickLabel', secondaryLabels,'HitTest','off',...
                'XTickLabelRotation',0,'box','off','FontSize',13,...%h.axis.FontSize,...
                'FontWeight',h.axis.FontWeight,'FontName',h.axis.FontName);
            h.axis.(posProp)([2,4]) = h.axis2.(posProp)([2,4]); % make room in original axes for 2ndary labels.
            h.axis2.(posProp)([1,3]) = h.axis.(posProp)([1,3]); % let original axis control lateral placement
            h.axis2.UserData.hlink = linkprop([h.axis, h.axis2],...
                {'Units',posProp,'ActivePositionProperty','Parent'}); % [5]
            linkaxes([h.axis, h.axis2], 'xy')
            if ~isUIFigure % [4]
                uistack(h.axis2, 'down')
            end
            if isprop(h.axis2, 'Toolbar')
                h.axis2.Toolbar.Visible = 'off'; % ver >= r2018b
            end
            h.axis2.XRuler.Axle.Visible = 'off';
            h.axis2.YRuler.Axle.Visible = 'off';
            h.axis.DeleteFcn = @(~,~)delete(h.axis2); % Delete axis2 if main axis is deleted
            set(h.figure,'CurrentAxes',h.axis)
        end
    end
    if any(strcmpi(groupLabelType, {'vertical','both'})) && ~verLessThan('Matlab','9.5') % r18b
        spaces = setdiff(1-p.Results.interGroupSpace : maxX, [x{:}]);
        endSpaceIdx = [diff(spaces),2] > 1;
        midSpace = spaces(endSpaceIdx) - (p.Results.interGroupSpace-1)/2;
        h.xline = arrayfun(@(x)xline(h.axis, x,'FontSize',h.axis.FontSize,...
            'FontWeight',h.axis.FontWeight,'FontName',h.axis.FontName),midSpace);
        set(h.xline(:), {'Label'}, cellstr(p.Results.secondarylabels(:))) % cellstr in case lbls are str
        vertLinesDrawn = true;
    end
end
% Draw vertical lines if requested and if they don't already exist.
if p.Results.groupLines && ~vertLinesDrawn && ~verLessThan('Matlab','9.5') %r18b
    spaces = setdiff(1:maxX+p.Results.interGroupSpace, [x{:}]);
    endSpaceIdx = [diff(spaces),2] > 1;
    midSpace = spaces(endSpaceIdx) - (p.Results.interGroupSpace-1)/2;
    h.xline = arrayfun(@(x)xline(h.axis, x,'-k'),midSpace);
end
clear('returnHoldState','TLOcleanup')
%% Return output only if requested
if nargout>0
    % Get and organize new boxplot groups
    bpobjPost = findobj(h.axis,'tag',bptag);
    h.boxplotGroup = bpobjPost(~ismember(bpobjPost, bpobjPre));
    handles = h;
end
end
function c = colorexpansion(colors, p, nGroups, nMembers)
% colors is color data. As of r2021a, boxplot 'Colors' can be RGB triplet/matrix
%   char vec, or string scalar of chars ("rgb"). Long color names is not accepted
%   by boxplot. 'colors' cannot be empty for this function.
% c: if 'colors' specifies more than 1 color, c is the color scheme expanded according
%   to GroupType. Otherwise, c is the same as colors.
% Other inputs defined in main func.
if isnumeric(colors) && size(colors,1)>1
    basecolors = colors;

elseif (ischar(colors) || isa(colors,'string')) && numel(char(colors))>1
    basecolors = char(colors);
    basecolors = basecolors(:); % col vec

else
    % If colors is not numeric, char, or string let boxplot throw the error.
    % If colors specifies only 1 color, copy colors to output.
    c = colors;
    return
end
isBetweenGroups = strcmpi('betweenGroups', p.Results.GroupType);
n = size(basecolors,1);
getRowIdx = @(n,m)round(mod(1:n,m+1E-08));
if isBetweenGroups
    % The first nMembers of colors will be used
    % Let boxplot do the expansion.
    rowNum = getRowIdx(nMembers,n);
    c = [basecolors(rowNum,:);repmat(basecolors(1,:),p.Results.interGroupSpace,1)];
else
    % The first nGroups colors will be used
    rowNum = getRowIdx(nGroups,n);
    c = repelem(basecolors(rowNum,:),nMembers+p.Results.interGroupSpace,1);
end
if ischar(c)
    c = c';
end
end
