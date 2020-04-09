require 'fileutils'

L = "/homes/users/asantini/local/lib:/homes/users/asantini/local/lib64"
G = "/homes/users/asantini/.gurobi/$HOSTNAME/gurobi.lic"
E = "/homes/users/asantini/local/src/tbkp/build/tbkp"
S = <<~EOF
    #!/bin/bash
    #SBATCH --partition=normal
    #SBATCH --time=02:00:00
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=1
    #SBATCH --cpus-per-task=1
    #SBATCH --mem-per-cpu=4GB
EOF

def create_script(instance, use_de, use_boole)
    instance = File.join(
        '/homes/users/asantini/local/src/tbkp/data/generated-instances',
        File.basename(instance)
    )

    b = File.basename(instance, '.txt')
    b += "-d" if use_de
    b += "-b" if use_boole

    script_f = File.join('scripts', "launch-#{b}.sh")
    error_f = File.join('/homes/users/asantini/local/src/tbkp/cluster/scripts', "err-#{b}.txt")
    output_f = File.join('/homes/users/asantini/local/src/tbkp/cluster/scripts', "out-#{b}.txt")
    results_f = File.join('/homes/users/asantini/local/src/tbkp/cluster/output', "res-#{b}.txt")
    
    params = ""
    params += " -d" if use_de
    params += " -b" if use_boole

    script = <<~EOF
        #{S.strip}
        #SBATCH -o #{output_f}
        #SBATCH -e #{error_f}

        module load Gurobi/9.0.0-lic
        LD_LIBRARY_PATH=#{L} GRB_LICENSE_FILE=#{G} #{E} -i #{instance} -o #{results_f} -t 3600#{params}
    EOF

    File.write(script_f, script)
end

def create_all_scripts
    Dir.glob('../data/generated-instances/*.txt') do |instance|
        [true, false].each do |use_de|
            [true, false].each do |use_boole|
                create_script(instance, use_de, use_boole)
            end
        end
    end
end

FileUtils.mkdir_p('scripts')
FileUtils.mkdir_p('output')
create_all_scripts