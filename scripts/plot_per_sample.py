import os
import subprocess

def get_sample_ids():
    fp = '/scratch/chd5n/aneuploidy/results/tables/'
    fn = 'samples_sd_rank.tsv'
    id_list = []
    with open(fp+fn,'r') as infile:
        line_num = 0
        for inline in infile:
            line_num = line_num + 1
            if line_num > 1:
                id_list.append(inline.strip('\n').split('\t')[0])
    return id_list

def shell_cmd(script_home,id_list):
    os.chdir(script_home)
    for each in id_list:
        cmd = ' '.join(['python run_plotter.py --subject_id', each])
        subprocess.call(cmd, shell=True)

def main():
    home = os.environ['HOME']
    script_home = home + '/projects/aneuploidy/scripts/'
    id_list = get_sample_ids()
    shell_cmd(script_home,id_list)

if __name__ == '__main__':
    main()
