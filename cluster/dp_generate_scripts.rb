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
    #SBATCH --mem-per-cpu=8GB
EOF

def create_script(instance, timeout: 3600)
    instance = File.join(
        '/homes/users/asantini/local/src/tbkp/data/generated-instances',
        File.basename(instance)
    )

    b = File.basename(instance, '.txt')
    script_f = File.join('scripts', "launch-#{b}.sh")
    error_f = File.join('/homes/users/asantini/local/src/tbkp/cluster/scripts', "err-#{b}.txt")
    output_f = File.join('/homes/users/asantini/local/src/tbkp/cluster/scripts', "out-#{b}.txt")
    results_f = File.join('/homes/users/asantini/local/src/tbkp/cluster/output', "res-#{b}.txt")

    script = <<~EOF
        #{S.strip}
        #SBATCH -o #{output_f}
        #SBATCH -e #{error_f}

        module load Gurobi/9.0.0-lic
        LD_LIBRARY_PATH=#{L} #{E} -s dp -i #{instance} -o #{results_f} -t #{timeout}
    EOF

    File.write(script_f, script)
end

def create_scripts(type, size)
    Dir.glob("../data/generated-instances/#{type}-#{size}-*.txt") do |instance|
        create_script(instance)
    end
end

unless ARGV.size == 2
    puts "Usage: dp_generate_scripts.rb [INSTANCE_TYPE] [INSTANCE_SIZE]"
end

unless ['type1', 'type2', 'type3', 'type4', 'type5'].include? ARGV[0]
    puts "Wrong instance type: #{ARGV[0]}"
end

unless [100, 500, 1000, 5000].include? ARGV[1].to_i
    puts "Wrong instance size: #{ARGV[1]}"
end

FileUtils.mkdir_p('scripts')
FileUtils.mkdir_p('output')

create_scripts(ARGV[0], ARGV[1])