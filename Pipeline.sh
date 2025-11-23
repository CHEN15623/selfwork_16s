#!/bin/bash
# ------------------------------------------------------------------
# 试运行环境windows11,r-3.4.5
# ------------------------------------------------------------------

# === 1. 全局配置 (可修改此处) ===
set -e                          # 遇到错误立即停止
WD="/d/16s"                     # 工作目录
DB="/d/16s/db"                  # 数据库目录
THREADS=4                       # 线程数
GROUP_COL="Group"               # Metadata中的分组列名
COMPARE_PAIR="KO-WT"            # 差异比较组 (格式: 实验-对照)
SEED=1                          # 随机种子

# 环境初始化
export PATH=$PATH:${DB}/win
cd "${WD}" || { echo "目录不存在: ${WD}"; exit 1; }

# 创建所有必要的输出目录
mkdir -p temp \
    result/raw result/alpha result/beta result/tax \
    result/gg result/compare result/stamp result/lefse

echo ">>> [1/10] 数据准备与环境检查..."
# 格式化元数据
sed -i 's/\r//' result/metadata.txt
# 检查并解压数据库(若不存在)
for db_file in "${DB}/usearch/rdp_16s_v18.fa" "${DB}/gg/97_otus.fa"; do
    if [ ! -f "${db_file}" ]; then
        gunzip -c "${db_file}.gz" > "${db_file}"
    fi
done

echo ">>> [2/10] 序列处理 (Merge & Rename)..."
# 使用 rush 并行合并双端序列
tail -n+2 result/metadata.txt | cut -f 1 | \
rush -j ${THREADS} "vsearch --fastq_mergepairs seq/{}_1.fq.gz --reverse seq/{}_2.fq.gz \
    --fastqout temp/{}.merged.fq --relabel {}."
# 合并为单个文件
cat temp/*.merged.fq > temp/all.fq

echo ">>> [3/10] 质控与聚类 (QC & Clustering)..."
# 切除引物与质控
vsearch --fastx_filter temp/all.fq \
    --fastq_stripleft 29 --fastq_stripright 18 \
    --fastq_maxee_rate 0.01 --fastaout temp/filtered.fa
# 去冗余
vsearch --derep_fulllength temp/filtered.fa \
    --minuniquesize 10 --sizeout --relabel Uni_ --output temp/uniques.fa
# Denoise (ASV)
usearch -unoise3 temp/uniques.fa -minsize 10 -zotus temp/zotus.fa
sed 's/Zotu/ASV_/g' temp/zotus.fa > temp/otus.fa
# (可选) 去嵌合体
vsearch --uchime_ref temp/otus.fa -db ${DB}/usearch/rdp_16s_v18.fa --nonchimeras result/raw/otus.fa
sed -i 's/\r//g' result/raw/otus.fa

echo ">>> [4/10] 生成特征表与过滤 (Feature Table)..."
# 映射生成 OTU 表
vsearch --usearch_global temp/filtered.fa -db result/raw/otus.fa \
    --id 0.97 --threads ${THREADS} --otutabout result/raw/otutab.txt
sed -i 's/\r//' result/raw/otutab.txt

# 物种注释
vsearch --sintax result/raw/otus.fa -db ${DB}/usearch/rdp_16s_v18.fa \
    --sintax_cutoff 0.1 --tabbedout result/raw/otus.sintax
sed -i 's/\r//' result/raw/otus.sintax

# 过滤非细菌/线粒体/叶绿体
Rscript ${DB}/script/otutab_filter_nonBac.R \
    --input result/raw/otutab.txt --taxonomy result/raw/otus.sintax \
    --output result/otutab.txt --stat result/raw/otutab_nonBac.stat --discard result/raw/otus.sintax.discard

# 更新 Fasta 和 Annotation
cut -f 1 result/otutab.txt | tail -n+2 > result/otutab.id
usearch -fastx_getseqs result/raw/otus.fa -labels result/otutab.id -fastaout result/otus.fa
awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}' result/raw/otus.sintax result/otutab.id > result/otus.sintax

# 等量抽平
Rscript ${DB}/script/otutab_rare.R --input result/otutab.txt \
    --depth 10000 --seed ${SEED} --normalize result/otutab_rare.txt --output result/alpha/vegan.txt

echo ">>> [5/10] Alpha 多样性分析与绘图..."
# 计算指标
usearch -alpha_div result/otutab_rare.txt -output result/alpha/alpha.txt
usearch -alpha_div_rare result/otutab_rare.txt -output result/alpha/alpha_rare.txt -method without_replacement
sed -i "s/-/\t0.0/g" result/alpha/alpha_rare.txt

# 1. 箱线图 (循环绘制所有指数)
for idx in $(head -n1 result/alpha/vegan.txt | cut -f 2-); do
    Rscript ${DB}/script/alpha_boxplot.R --alpha_index ${idx} \
    --input result/alpha/vegan.txt --design result/metadata.txt \
    --group ${GROUP_COL} --output result/alpha/ --width 89 --height 59
done
# 2. 柱状图
Rscript ${DB}/script/alpha_barplot.R --alpha_index richness \
    --input result/alpha/vegan.txt --design result/metadata.txt \
    --group ${GROUP_COL} --output result/alpha/ --width 89 --height 59
# 3. 稀释曲线
Rscript ${DB}/script/alpha_rare_curve.R \
    --input result/alpha/alpha_rare.txt --design result/metadata.txt \
    --group ${GROUP_COL} --output result/alpha/ --width 120 --height 59
# 4. 维恩图 (需生成按组均值文件)
Rscript ${DB}/script/otu_mean.R --input result/otutab.txt --metadata result/metadata.txt \
    --group ${GROUP_COL} --thre 0 --scale TRUE --zoom 100 --all TRUE --type mean --output result/otutab_mean.txt
# 筛选>0.05的OTU用于Venn
awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i; print "OTU","Group";} \
    else {for(i=2;i<=NF;i++) if($i>0.05) print $1, a[i];}}' result/otutab_mean.txt > result/alpha/otu_group_exist.txt
# 绘制 Venn (根据你的分组修改 -a -b -c 参数)
bash ${DB}/script/sp_vennDiagram.sh -f result/alpha/otu_group_exist.txt \
    -a WT -b KO -c OE -w 3 -u 3 -p WT_KO_OE

echo ">>> [6/10] Beta 多样性分析与绘图..."
# 建树与计算距离矩阵
usearch -cluster_agg result/otus.fa -treeout result/otus.tree
usearch -beta_div result/otutab_rare.txt -tree result/otus.tree -filename_prefix result/beta/

# 1. PCoA (无标签 & 有标签)
for label_flag in FALSE TRUE; do
    out_suffix=$( [ "$label_flag" == "TRUE" ] && echo ".label" || echo "" )
    Rscript ${DB}/script/beta_pcoa.R --input result/beta/bray_curtis.txt --design result/metadata.txt \
        --group ${GROUP_COL} --label ${label_flag} --width 89 --height 59 \
        --output result/beta/bray_curtis.pcoa${out_suffix}.pdf
done

# 2. CPCoA (限制性主坐标)
Rscript ${DB}/script/beta_cpcoa.R --input result/beta/bray_curtis.txt --design result/metadata.txt \
    --group ${GROUP_COL} --label TRUE --width 89 --height 59 --output result/beta/bray_curtis.cpcoa.label.pdf

# 3. 热图 (Pheatmap)
cut -f 1-2 result/metadata.txt > temp/group.txt
bash ${DB}/script/sp_pheatmap.sh -f result/beta/bray_curtis.txt -H 'TRUE' -u 6.9 -v 5.6 \
    -P temp/group.txt -Q temp/group.txt

echo ">>> [7/10] 物种分类汇总与绘图 (Taxonomy)..."
# 格式化 Taxonomy 文件 (8列格式)
cut -f 1,4 result/otus.sintax | sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g' > result/taxonomy2.txt
awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
      split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
      print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' result/taxonomy2.txt > temp/otus.tax
sed 's/;/\t/g;s/.__//g;' temp/otus.tax | cut -f 1-8 | \
sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' > result/taxonomy.txt

# 循环处理 p,c,o,f,g 各级别
for rank in p c o f g; do
    # 汇总
    usearch -sintax_summary result/otus.sintax -otutabin result/otutab_rare.txt -rank ${rank} -output result/tax/sum_${rank}.txt
    sed -i 's/(//g;s/)//g;s/\"//g;s/\#//g;s/\/Chloroplast//g' result/tax/sum_${rank}.txt
    
    # 1. 堆叠柱状图 (Stackplot)
    Rscript ${DB}/script/tax_stackplot.R --input result/tax/sum_${rank}.txt --design result/metadata.txt \
        --group ${GROUP_COL} --output result/tax/sum_${rank}.stackplot --legend 8 --width 89 --height 59
    
    # 2. 弦图 (Circlize, 前5组)
    if [ "$rank" == "c" ]; then # 示例仅对Class做弦图，避免太多
        Rscript ${DB}/script/tax_circlize.R --input result/tax/sum_${rank}.txt --design result/metadata.txt --group ${GROUP_COL} --legend 5
        mv circlize.pdf result/tax/sum_${rank}.circlize.pdf
        mv circlize_legend.pdf result/tax/sum_${rank}.circlize_legend.pdf
    fi
done

# 3. 树图 (Treemap)
Rscript ${DB}/script/tax_maptree.R --input result/otutab.txt --taxonomy result/taxonomy.txt \
    --output result/tax/tax_maptree.pdf --topN 100 --width 183 --height 118

echo ">>> [8/10] GreenGenes 有参比对 (Functional Prediction Prep)..."
usearch -otutab temp/filtered.fa -otus ${DB}/gg/97_otus.fa \
    --otutabout result/gg/otutab.txt --threads ${THREADS}
usearch -otutab_stats result/gg/otutab.txt -output result/gg/otutab.stat

echo ">>> [9/10] 差异分析与 STAMP/LEfSe 准备..."
# 1. 差异比较与火山图
Rscript ${DB}/script/compare.R --input result/otutab.txt --design result/metadata.txt \
    --group ${GROUP_COL} --compare ${COMPARE_PAIR} --threshold 0.1 --method edgeR --output result/compare/
Rscript ${DB}/script/compare_volcano.R --input result/compare/${COMPARE_PAIR}.txt \
    --output result/compare/${COMPARE_PAIR}.volcano.pdf --width 89 --height 59

# 2. 差异热图
bash ${DB}/script/compare_heatmap.sh -i result/compare/${COMPARE_PAIR}.txt -l 7 \
    -d result/metadata.txt -A ${GROUP_COL} -t result/taxonomy.txt -w 8 -h 5 -s 7 -o result/compare/${COMPARE_PAIR}

# 3. 曼哈顿图 (针对 Phylum 和 Genus)
bash ${DB}/script/compare_manhattan.sh -i result/compare/${COMPARE_PAIR}.txt -t result/taxonomy.txt \
    -p result/tax/sum_p.txt -w 250 -v 59 -s 7 -l 10 -o result/compare/${COMPARE_PAIR}.manhattan.p.pdf
bash ${DB}/script/compare_manhattan.sh -i result/compare/${COMPARE_PAIR}.txt -t result/taxonomy.txt \
    -p result/tax/sum_g.txt -w 250 -v 59 -s 7 -l 10 -L Genus -o result/compare/${COMPARE_PAIR}.manhattan.g.pdf

# 4. STAMP 准备与扩展柱状图
Rscript ${DB}/script/format2stamp.R --input result/otutab.txt --taxonomy result/taxonomy.txt \
    --threshold 0.01 --output result/stamp/tax
# 绘制 STAMP 风格图 (属水平)
Rscript ${DB}/script/compare_stamp.R --input result/stamp/tax_5Family.txt --metadata result/metadata.txt \
    --group ${GROUP_COL} --compare ${COMPARE_PAIR} --threshold 0.001 --method "t.test" --width 280 --height 159 \
    --output result/stamp/${COMPARE_PAIR}

# 5. LEfSe 输入文件生成
Rscript ${DB}/script/format2lefse.R --input result/otutab.txt --taxonomy result/taxonomy.txt \
    --design result/metadata.txt --group ${GROUP_COL} --threshold 0.4 --output result/lefse/LEfSe

echo ">>> [10/10] 清理与 MD5 校验..."
rm -rf temp/*.fq temp/*.fa
cd seq
md5sum *_1.fq.gz > md5sum1.txt
md5sum *_2.fq.gz > md5sum2.txt
paste md5sum1.txt md5sum2.txt | awk '{print $2"\t"$1"\t"$4"\t"$3}' | sed 's/*//g' > ../result/md5sum.txt
rm md5sum*
cd ..

echo ">>> 分析全部完成！请查看 result 目录。"
