touch output.txt
python kosmatau3d/terminal_scripts/run_model_grid.py -o true -m 8 --grid cm-f_fuv |& tee -a output.txt
python kosmatau3d/terminal_scripts/run_model_grid.py -o true -m 8 --grid convergence |& tee  output.txt
python kosmatau3d/terminal_scripts/run_model_grid.py -o true -m 8 --grid f_cm-cm |& tee -a output.txt
python kosmatau3d/terminal_scripts/run_model_grid.py -o true -m 8 --grid f_cm-f_icm |& tee -a output.txt
python kosmatau3d/terminal_scripts/run_model_grid.py -o true -m 8 --grid f_hi-f_fuv |& tee -a output.txt
python kosmatau3d/terminal_scripts/run_model_grid.py -o true -m 8 --grid f_fuv_gc |& tee -a output.txt
python kosmatau3d/terminal_scripts/run_model_grid.py -o true -m 8 --grid fuv_cl |& tee -a output.txt
python kosmatau3d/terminal_scripts/run_model_grid.py -o true -m 8 --grid r_cmz-f_fuv |& tee -a output.txt
python kosmatau3d/terminal_scripts/run_model_grid.py -o true -m 8 --grid vox_disp |& tee -a output.txt
python kosmatau3d/terminal_scripts/run_model_grid.py -o true -m 8 --grid wnm |& tee -a output.txt
