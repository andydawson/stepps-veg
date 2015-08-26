
import subprocess

runs = [ ('veg_120knots_c1',
          './veg_od_mpp_262.exe \
          sample num_warmup=150 num_samples=1000 save_warmup=1\
          data file=../r/dump/veg_data_12taxa_6341cells_120knots_v0.4.dump \
          init=../r/dump/veg_data_12taxa_6341cells_120knots_v0.4_inits.dump \
          output file=../output/12taxa_6341cells_120knots_v0.4_c1.csv\
          random seed=42'),
         ('veg_120knots_c2',
          './veg_od_mpp_262.exe \
          sample num_warmup=150 num_samples=1000 save_warmup=1\
          data file=../r/dump/veg_data_12taxa_6341cells_120knots_v0.4.dump \
          init=../r/dump/veg_data_12taxa_6341cells_120knots_v0.4_inits.dump \
          output file=../output/12taxa_6341cells_120knots_v0.4_c2.csv\
          random seed=42')
]

qsub = """\
#!/bin/sh
#SBATCH --job-name={name}
#SBATCH -c {threads}
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andria.dawson@gmail.com
#SBATCH --exclude sm-20

cd $HOME/Documents/projects/stepps-veg/cpp
export OMP_NUM_THREADS={threads}
srun {command}
"""

dry_run = False

for name, command in runs:
#    sub = qsub.format(queue="low.q", walltime="672:00:00", command=run, threads=1)
    sub = qsub.format(command=command, threads=12, name=name)
    with open(name + ".sh", 'w') as f:
        f.write(sub)
    print "submitting:", name
    if not dry_run:
        subprocess.check_call(['sbatch', name + '.sh'])
    else:
        print sub
