require 'fileutils'

MEM = (ARGV[0] == 'mem' ? '8GB' : '4GB')

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
    #SBATCH --mem-per-cpu=#{MEM}
EOF

def create_script(instance, early_combo, use_de, use_boole, boole_freq)
    instance = File.join(
        '/homes/users/asantini/local/src/tbkp/data/generated-instances',
        File.basename(instance)
    )

    b = File.basename(instance, '.txt')
    b += "-d" if use_de
    b += "-b#{boole_freq}" if use_boole
    b += "-c" if early_combo

    script_f = File.join('scripts', "launch-#{b}.sh")
    error_f = File.join('/homes/users/asantini/local/src/tbkp/cluster/scripts', "err-#{b}.txt")
    output_f = File.join('/homes/users/asantini/local/src/tbkp/cluster/scripts', "out-#{b}.txt")
    results_f = File.join('/homes/users/asantini/local/src/tbkp/cluster/output', "res-#{b}.txt")
    
    params = ""
    params += " -d 1" if use_de
    params += " -b 1 -f #{boole_freq}" if use_boole
    params += " -c 1" if early_combo

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
                next if use_boole && !use_de # Doesn't make sense
                [true, false].each do |early_combo|
                    create_script(instance, early_combo, use_de, use_boole, 1)
                    create_script(instance, early_combo, use_de, use_boole, 100) if use_boole
                end
            end
        end
    end
end

FileUtils.mkdir_p('scripts')
FileUtils.mkdir_p('output')
create_all_scripts
