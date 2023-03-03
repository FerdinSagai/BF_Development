pwd=` pwd`
case_name="Case_1"
input_folder_name="Input_"$case_name
output_folder_name="Result_"$case_name
export OMP_NUM_THREADS=4

rm -rf $pwd/$output_folder_name
mkdir $pwd/$output_folder_name

cd $pwd/$input_folder_name
cp $pwd/$input_folder_name/1.*	    $pwd/$output_folder_name


cd $pwd/inputscript_reader
make

cd $pwd/src
make

mv $pwd/inputscript_reader/INPUT_GUI	$pwd/$output_folder_name
mv $pwd/src/GF_SOLVER					$pwd/$output_folder_name

cd $pwd/$output_folder_name
./INPUT_GUI
./GF_SOLVER
#nohup ./GF_SOLVER &