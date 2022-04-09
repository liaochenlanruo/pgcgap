#!/bin/sh
# mkdir -p $PREFIX/bin
reMatch() {
	typeset ec
	unset -v reMatch # initialize output variable
	[[ $1 =~ $2 ]] # perform the regex test
	ec=$? # save exit code
	if [[ $ec -eq 0 ]]; then # copy result to output variable
		[[ -n $BASH_VERSION ]] && reMatch=( "${BASH_REMATCH[@]}" )
		[[ -n $KSH_VERSION ]]  && reMatch=( "${.sh.match[@]}" )
		[[ -n $ZSH_VERSION ]]  && reMatch=( "$MATCH" "${match[@]}" )
	fi
	return $ec
}

dir=`which prokka`

reMatch $dir "(^\S+)\/bin\/prokka"
PREFIX=${reMatch[1]}

wget https://codeload.github.com/liaochenlanruo/pgcgap/tar.gz/refs/tags/v1.0.34
tar zxvf v1.0.34
cd pgcgap-1.0.34

cp pgcgap.pl $PREFIX/bin/pgcgap

cp Functions/Pan/plot_3Dpie.R $PREFIX/bin/
cp Functions/Pan/fmplot.py $PREFIX/bin/
cp Functions/Pan/grep_cds_aas_from_gff3.pl $PREFIX/bin/

cp Functions/COG/COG.pl $PREFIX/bin/
cp Functions/COG/COG2020.pl $PREFIX/bin/
cp Functions/COG/COGdiamond2022.pl $PREFIX/bin/
cp Functions/COG/get_flag_relative_abundances_table.pl $PREFIX/bin/
cp Functions/COG/Plot_COG.R $PREFIX/bin/
cp Functions/COG/Plot_COG_Abundance.R $PREFIX/bin/

cp Functions/ANI/triangle2list.pl $PREFIX/bin/
cp Functions/ANI/get_ANImatrix.pl $PREFIX/bin/
cp Functions/ANI/Plot_ANIheatmap.R $PREFIX/bin/

cp Functions/MASH/get_Mash_Matrix.pl $PREFIX/bin/
cp Functions/MASH/Plot_MashHeatmap.R $PREFIX/bin/

cp Functions/Assemble/genome_LenFilter_stats.pl $PREFIX/bin/
cp Functions/Assemble/get_stats_summary.pl $PREFIX/bin/

chmod a+x $PREFIX/bin/pgcgap
chmod a+x $PREFIX/bin/plot_3Dpie.R
chmod a+x $PREFIX/bin/fmplot.py
chmod a+x $PREFIX/bin/grep_cds_aas_from_gff3.pl
chmod a+x $PREFIX/bin/COG.pl
chmod a+x $PREFIX/bin/COG2020.pl
chmod a+x $PREFIX/bin/COGdiamond2022.pl
chmod a+x $PREFIX/bin/get_flag_relative_abundances_table.pl
chmod a+x $PREFIX/bin/Plot_COG.R
chmod a+x $PREFIX/bin/triangle2list.pl
chmod a+x $PREFIX/bin/get_ANImatrix.pl
chmod a+x $PREFIX/bin/Plot_ANIheatmap.R
chmod a+x $PREFIX/bin/Plot_COG_Abundance.R
chmod a+x $PREFIX/bin/get_Mash_Matrix.pl
chmod a+x $PREFIX/bin/Plot_MashHeatmap.R
chmod a+x $PREFIX/bin/genome_LenFilter_stats.pl
chmod a+x $PREFIX/bin/get_stats_summary.pl

pgcgap --version
