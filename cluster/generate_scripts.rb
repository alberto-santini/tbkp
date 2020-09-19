require 'fileutils'

MEM = (ARGV[0] == 'mem' ? '16GB' : '8GB')

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

def create_script(instance, early_combo:, use_de:, use_boole:, boole_freq:, use_cr:, all_bounds: false, max_nodes: 0)
    instance = File.join(
        '/homes/users/asantini/local/src/tbkp/data/generated-instances',
        File.basename(instance)
    )

    b = File.basename(instance, '.txt')
    b += "-d" if use_de
    b += "-b#{boole_freq}" if use_boole
    b += "-r" if use_cr
    b += "-c" if early_combo

    script_f = File.join('scripts', "launch-#{b}.sh")
    error_f = File.join('/homes/users/asantini/local/src/tbkp/cluster/scripts', "err-#{b}.txt")
    output_f = File.join('/homes/users/asantini/local/src/tbkp/cluster/scripts', "out-#{b}.txt")
    results_f = File.join('/homes/users/asantini/local/src/tbkp/cluster/output', "res-#{b}.txt")
    
    params = ""
    params += " -d 1" if use_de
    params += " -b 1 -f #{boole_freq}" if use_boole
    params += " -c 1" if early_combo
    params += " -r 1" if use_cr
    params += " -a 1" if all_bounds
    params += " -n #{max_nodes}"

    script = <<~EOF
        #{S.strip}
        #SBATCH -o #{output_f}
        #SBATCH -e #{error_f}

        module load Gurobi/9.0.0-lic
        LD_LIBRARY_PATH=#{L} #{E} -i #{instance} -o #{results_f} -t 3600#{params}
    EOF

    File.write(script_f, script)
end

def create_boole_root_bound_scripts
    Dir.glob('../data/generated-instances/*.txt') do |instance|
        create_script instance, early_combo: true, use_de: false, use_boole: true, boole_freq: 1, use_cr: false, all_bounds: false, max_nodes: 1
    end
end

def create_bb_eval_scripts
    Dir.glob('../data/generated-instances/*.txt') do |instance|
        # First configuration: enumeration
        # create_script instance, early_combo: true, use_de: false, use_boole: false, boole_freq: 0, use_cr: false
        # Second configuration: DEbounds (z1lower + d1upper)
        # create_script instance, early_combo: true, use_de: true, use_boole: false, boole_freq: 0, use_cr: false
        # Third configuration: all bounds, i.e., DEbounds (z1lower + z1upper), BOOLEbound (z2lower), CRbound (z2upper)
        # create_script instance, early_combo: true, use_de: true, use_boole: true, boole_freq: 1, use_cr: true
        # Fourth configuration: DEbounds (z1lower + d1upper) + BOOLEbound (z2lower)
        # create_script instance, early_combo: true, use_de: true, use_boole: true, boole_freq: 1, use_cr: false
        # Fifth configuration: DEbounds (z1lower + d1upper) + CRbound (z2upper)
        # create_script instance, early_combo: true, use_de: true, use_boole: false, boole_freq: 0, use_cr: true
        # Sixth configuration: DEbounds (z1lower + d1upper) + CRbound (z2upper) but *without* early combo
        create_script instance, early_combo: false, use_de: true, use_boole: false, boole_freq: 0, use_cr: true
    end
end

FileUtils.mkdir_p('scripts')
FileUtils.mkdir_p('output')

create_bb_eval_scripts