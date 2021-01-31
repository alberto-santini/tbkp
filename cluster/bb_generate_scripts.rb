require 'fileutils'

L = "/homes/users/asantini/local/lib:/homes/users/asantini/local/lib64"
E = "/homes/users/asantini/local/src/tbkp/build/tbkp"
S = <<~EOF
    #!/bin/bash
    #SBATCH --partition=normal
    #SBATCH --time=02:00:00
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=1
    #SBATCH --cpus-per-task=1
EOF

def create_script(
    instance,
    config_n,
    early_combo: true,
    early_pruning: true,
    use_de: false,
    use_boole: false,
    use_cr: false,
    boole_freq: 1,
    boole_tl: 1,
    boole_root_tl: 1,
    max_nodes: 0
)
    instance = File.join(
        '/homes/users/asantini/local/src/tbkp/data/generated-instances',
        File.basename(instance)
    )

    f = File.basename(instance, '.txt')
    sz = f.split('-')[1].to_i

    b = "config-#{config_n}-#{f}"
    b += "-d" if use_de
    b += "-b#{boole_freq},#{boole_tl},#{boole_root_tl}" if use_boole
    b += "-r" if use_cr
    b += "-c" if early_combo
    b += "-p" if early_pruning

    results_base = '/homes/users/asantini/local/src/tbkp/cluster/output'
    FileUtils.mkdir_p(File.join(results_base, "bb-config-#{config_n}"))

    script_f = File.join('scripts', "launch-bb-#{b}.sh")
    error_f = File.join('/homes/users/asantini/local/src/tbkp/cluster/scripts', "err-#{b}.txt")
    output_f = File.join('/homes/users/asantini/local/src/tbkp/cluster/scripts', "out-#{b}.txt")
    results_f = File.join(results_base, "bb-config-#{config_n}", "res-#{b}.txt")
    
    params = ""
    params += (use_de ? " -d 1" : " -d 0")
    params += (use_boole ? " -b 1 -f #{boole_freq} -T #{boole_tl} -R #{boole_root_tl}" : " -b 0")
    params += (use_cr ? " -r 1" : " -r 0")
    params += (early_combo ? " -c 1" : " -c 0")
    params += (early_pruning ? " -p 1" : " -p 0")
    params += " -n #{max_nodes}"

    mem = case sz
        when 5000 then '32GB'
        when 1000 then '8GB'
        else '4GB'
    end

    script = <<~EOF
        #{S.strip}
        #SBATCH --mem-per-cpu=#{mem}
        #SBATCH -o #{output_f}
        #SBATCH -e #{error_f}

        module load Gurobi/9.0.0-lic
        LD_LIBRARY_PATH=#{L} #{E} -i #{instance} -o #{results_f} -t 3600#{params}
    EOF

    File.write(script_f, script)
end

def create_bb_eval_scripts(configuration)
    Dir.glob('../data/generated-instances/*.txt') do |instance|
        case configuration
        when 1
            # Enumeration
            create_script instance, 1, early_pruning: false
        when 2
            # DEbounds (z1lower + z1upper)
            create_script instance, 2, use_de: true
        when 3
            # DEbounds (z1lower + z1upper) + CRbound (z2upper)
            create_script instance, 3, use_de: true, use_cr: true
        when 4
            # DEbounds (z1lower + z1upper) + BOOLEbound (z2lower: 10s root, 1s other, at each node)
            create_script instance, 4, use_de: true,
                use_boole: true, boole_freq: 1, boole_tl: 1, boole_root_tl: 10
        when 5
            # DEbounds (z1lower + z1upper) + BOOLEbound (z2lower: 10s root, 1s other, each 100th node)
            create_script instance, 5, use_de: true,
                use_boole: true, boole_freq: 100, boole_tl: 1, boole_root_tl: 10
        when 6
            # DEbounds (z1lower + z1upper) + BOOLEbound (z2lower: 10s root, 1s other, each 1000th node)
            create_script instance, 5, use_de: true,
                use_boole: true, boole_freq: 1000, boole_tl: 1, boole_root_tl: 10
        when 7
            # DEbounds (z1lower + z1upper) + BOOLEbound (z2lower: 10s root, 1s other, each 100th node) + CRbound (z2upper)
            create_script instance, 6, use_de: true,
                use_boole: true, boole_freq: 100, boole_tl: 1, boole_root_tl: 10,
                use_cr: true
        when 8
            # DEbounds (z1lower + z1upper) + BOOLEbound (z2lower: 10s root, 1s other, each 1000th node) + CRbound (z2upper)
            create_script instance, 6, use_de: true,
                use_boole: true, boole_freq: 1000, boole_tl: 1, boole_root_tl: 10,
                use_cr: true
        when 9
            # DEbounds (z1lower + z1upper) + BOOLEbound (z2lower: 10s root, 1s other, each 1000th node) + CRbound (z2upper), no early combo
            create_script instance, 6, use_de: true,
                use_boole: true, boole_freq: 1000, boole_tl: 1, boole_root_tl: 10,
                use_cr: true, early_combo: false
        when 10
            # DEbounds (z1lower + z1upper) + BOOLEbound (z2lower: 10s root, 1s other, each 1000th node) + CRbound (z2upper), no early pruning
            create_script instance, 6, use_de: true,
                use_boole: true, boole_freq: 1000, boole_tl: 1, boole_root_tl: 10,
                use_cr: true, early_pruning: false
        end
    end
end

def create_bb_root_node_scripts(configuration)
    Dir.glob('../data/generated-instances/*.txt') do |instance|
        case configuration
        when 1
            # DEbounds (z1lower + z1upper)
            create_script instance, early_combo: false, early_pruning: false, use_de: true, use_boole: false, boole_freq: 0, boole_tl: 3600, quad_boole: true, use_cr: false, max_nodes: 1
        when 2
            # BOOLEbound (z2lower)
            create_script instance, early_combo: false, early_pruning: false, use_de: false, use_boole: true, boole_freq: 1, boole_tl: 3600, quad_boole: true, use_cr: false, max_nodes: 1
        when 3
            # BOOLEbound 1 second (z2lower)
            create_script instance, early_combo: false, early_pruning: false, use_de: false, use_boole: true, boole_freq: 1, boole_tl: 1, quad_boole: true, use_cr: false, max_nodes: 1
        when 4
            # BOOLEbound linearised 1 second (z2lower)
            create_script instance, early_combo: false, early_pruning: false, use_de: false, use_boole: true, boole_freq: 1, boole_tl: 1, quad_boole: false, use_cr: false, max_nodes: 1
        when 5
            # BOOLEbound 10 minutes (z2lower)
            create_script instance, early_combo: false, early_pruning: false, use_de: false, use_boole: true, boole_freq: 1, boole_tl: 600, quad_boole: true, use_cr: false, max_nodes: 1
        when 6
            # CRbound (z2upper)
            create_script instance, early_combo: false, early_pruning: false, use_de: false, use_boole: false, boole_freq: 0, boole_tl: 3600, quad_boole: true, use_cr: true, max_nodes: 1
        when 7
            # BOOLEbound 10 seconds (z2lower)
            create_script instance, early_combo: false, early_pruning: false, use_de: false, use_boole: true, boole_freq: 1, boole_tl: 10, quad_boole: true, use_cr: false, max_nodes: 1
        when 8
            # BOOLEbound 1 minute (z2lower)
            create_script instance, early_combo: false, early_pruning: false, use_de: false, use_boole: true, boole_freq: 1, boole_tl: 60, quad_boole: true, use_cr: false, max_nodes: 1
        end
    end
end

FileUtils.mkdir_p('scripts')
FileUtils.mkdir_p('output')

if ARGV[0] == 'root'
    create_bb_root_node_scripts(ARGV[1].to_i)
elsif ARGV[0] == 'bb'
    create_bb_eval_scripts(ARGV[1].to_i)
end
