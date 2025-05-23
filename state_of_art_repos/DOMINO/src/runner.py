import argparse
import os

import src.constants as constants
from src.core.domino import main as domino_main
from src.core.preprocess_slices import create_slices
from src.utils.visualize_modules import visualize_modules


def main_domino():

    parser = argparse.ArgumentParser(description='DOMINO: An active module identification algorithm with reduce rate of false.\n NOTE YOU SHOULD RUN THE SLICES SCRIPT FIRST! (more info, type slicer -h) \n Example input files are available @ https://github.com/Shamir-Lab/DOMINO/tree/master/examples')
    parser.add_argument('-a', '--active_genes_files', dest='active_genes_files', help='Comma delimited list of absolute paths to files, each containing a list of active genes, separated by a new line char (\\n). e.g. /path/to/active_genes_files_1,/path/to/active_genes_files_2.', default="examples/tnfa_active_genes_file.txt")
    parser.add_argument('-n', '--network_file', dest='network_file', help='A path to network file (sif format). e.g. /path/to/network_file.sif', default="examples/huri.sif")
    parser.add_argument('-s', '--slices_file', dest='slices_file', help='A path to slices file (i.e. the output of "slicer" script). e.g., /path/to/slices_file.txt', default="examples/huri_slices.txt")
    parser.add_argument('-o', '--output_folder', dest='output_folder', help='A folder where output files will be written e.g., /path/to/output', default="examples/output")
    parser.add_argument('-c', '--use_cache', dest='use_cache', help='Use auto-generated cache network files (*.pkl) from previous executions with the same network. NOTE: (1) THIS IS NOT THE SLICES FILE! (2) If the content of the file has changed, you should set this option to "false"', default="true")
    parser.add_argument('-p', '--parallelization', dest='parallelization', help='The number of threads allocated to the run (usually single thread is enough)', default="1")
    parser.add_argument('-v', '--visualization', dest='visualization', help='Indicates whether a visualization of the modules ought to be generated', default="true")
    parser.add_argument('-sth', '--slice_threshold', dest='slice_threshold', default="0.3", help='The threshold for considering a slice as relevant')
    parser.add_argument('-mth', '--module_threshold', dest='module_threshold', default="0.05", help='The threshold for considering a putative module as final module')


    args = parser.parse_args()
    active_genes_files = args.active_genes_files.split(",")
    output_folder = args.output_folder
    network_file = args.network_file
    slices_file = args.slices_file
    slice_threshold = float(args.slice_threshold)
    module_threshold = float(args.module_threshold)
    use_cache = args.use_cache=="true"
    parallelization = int(args.parallelization)
    visualization = args.visualization=="true"

    constants.N_OF_THREADS=parallelization
    constants.USE_CACHE=use_cache

    for cur_ag in active_genes_files:
        G_final_modules=domino_main(active_genes_file=cur_ag, network_file=network_file, slices_file=slices_file, slice_threshold=slice_threshold, module_threshold=module_threshold)
        activity_name=os.path.splitext(os.path.split(cur_ag)[-1])[0]
        report_folder=os.path.join(output_folder,activity_name)
        try:
            os.makedirs(report_folder)
        except:
            pass

        out_file=os.path.join(report_folder, "modules.out") 
        if len(G_final_modules) !=0:
            open(out_file, 'w+').write("\n".join(['[%s]' % ', '.join(list(m.nodes)) for m in G_final_modules])+"\n")
        else:
            open(out_file, 'w+').write("")

        print(f'{len(G_final_modules)} final modules are reported at {out_file}')
        print(visualization)
        if visualization:
            visualize_modules(os.path.splitext(cur_ag.split('/')[-1])[0], G_final_modules, None, network_file, report_folder)
        return G_final_modules

def main_slicer():

    parser = argparse.ArgumentParser(description='Slicer for DOMINO (step #0): A preprocessing step for the network')
    parser.add_argument('-n', '--network_file', dest='network_file', help='A path to network file (sif format). e.g. /path/to/network_file.sif', default="examples/huri.sif")
    parser.add_argument('-o', '--output_file', dest='output_file', default="examples/huri.sif", help='A path to the output slices file. e.g., /path/to/output/slices_file.txt')


    args = parser.parse_args()
    network_file = args.network_file
    output_file = args.output_file
    create_slices(network_file, output_file)




if __name__=="__main__":
    main_slicer()
    main_domino()
