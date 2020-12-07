require 'simple-random'

def generate_type1(n, b)
    f = Dir.glob(File.join('from-pisinger-hard', "knapPI_*_#{n}_*.csv")).to_a.sample
    i = rand(1..100)
    nn, capacity, weights, profits = get_pis_instance(f, i)

    raise "Wrong n from Pisinger instance (exp #{n}, got #{nn})!" if n != nn

    profits, weights = profits.zip(weights).sort.reverse.transpose

    n_timebombs = (n * b).ceil.to_i
    max_tb_profit = profits.max
    min_tb_profit = profits.take(n_timebombs).min
    max_det_profit = profits.select{|p| p < min_tb_profit}.max

    puts 'Warning: flat instance' if max_tb_profit == max_det_profit

    probabilities = Array.new(n_timebombs) do |i|
        1.0 - 0.1 * (profits[i] - max_det_profit) / (max_tb_profit - max_det_profit)
    end.map{|pi| pi.round(6)}

    (n - n_timebombs).times do
        probabilities << 1
    end

    return weights, profits, probabilities, capacity
end

def generate_type2(n, b)
    f = Dir.glob(File.join('from-pisinger-hard', "knapPI_*_#{n}_*.csv")).to_a.sample
    i = rand(1..100)
    nn, capacity, weights, profits = get_pis_instance(f, i)

    raise "Wrong n from Pisinger instance (exp #{n}, got #{nn})!" if n != nn

    densities = Array.new(n) do |i|
        profits[i].to_f / weights[i].to_f
    end
    densities, profits, weights = densities.zip(profits, weights).sort.reverse.transpose

    n_timebombs = (n * b).ceil.to_i
    max_tb_density = densities.max
    min_tb_density = densities.take(n_timebombs).min
    max_det_density = densities.select{|p| p < min_tb_density}.max || densities[n_timebombs]

    puts 'Warning: flat instance' if max_tb_density == max_det_density

    probabilities = Array.new(n_timebombs) do |i|
        1.0 - 0.1 * (densities[i] - max_det_density) / (max_tb_density - max_det_density)
    end.map{|pi| pi.round(6)}

    (n - n_timebombs).times do
        probabilities << 1
    end

    return weights, profits, probabilities, capacity
end

def generate_type3(n, b)
    f = Dir.glob(File.join('from-pisinger-hard', "knapPI_*_#{n}_*.csv")).to_a.sample
    i = rand(1..100)
    nn, capacity, weights, profits = get_pis_instance(f, i)

    raise "Wrong n from Pisinger instance (exp #{n}, got #{nn})!" if n != nn

    profits, weights = profits.zip(weights).sort.reverse.transpose

    r = SimpleRandom.new
    r.set_seed
    n_timebombs = (n * b).ceil.to_i
    probabilities = Array.new(n_timebombs) do
        1.0 - r.beta(1, 10)
    end.map{|pi| pi.round(6)}

    (n - n_timebombs).times do
        probabilities << 1
    end

    return weights, profits, probabilities, capacity
end

def generate_type4(n, b)
    f = Dir.glob(File.join('from-pisinger-hard', "knapPI_*_#{n}_*.csv")).to_a.sample
    i = rand(1..100)
    nn, capacity, weights, profits = get_pis_instance(f, i)

    raise "Wrong n from Pisinger instance (exp #{n}, got #{nn})!" if n != nn

    densities = Array.new(n) do |i|
        profits[i].to_f / weights[i].to_f
    end
    densities, profits, weights = densities.zip(profits, weights).sort.reverse.transpose

    r = SimpleRandom.new
    r.set_seed
    n_timebombs = (n * b).ceil.to_i
    probabilities = Array.new(n_timebombs) do
        1.0 - r.beta(1, 10)
    end.map{|pi| pi.round(6)}

    (n - n_timebombs).times do
        probabilities << 1
    end

    return weights, profits, probabilities, capacity
end

def generate_type5(n, k, c)
    w = (k * c / n).floor.to_i
    weights = Array.new(n) {w}

    r = SimpleRandom.new
    r.set_seed
    profits = Array.new(n) {[1, (w + r.normal(0, w/8.0)).ceil].max}.map(&:to_i)
    minprofit = profits.min

    probabilities = Array.new(n) do |i|
        minprofit.to_f / profits[i].to_f
    end.map{|pi| pi.round(6)}

    return weights, profits, probabilities
end

def get_pis_instance(f, i)
    ln = 0
    inst_n = 0

    n = 0
    c = 0
    weights = Array.new
    profits = Array.new

    File.readlines(f).each do |line|
        inst_n += 1 if (line.size > 5) && (line[0..5] == 'knapPI')
        
        if inst_n == i
            if ln == 1
                n = line.split[1].to_i
            elsif ln == 2
                c = line.split[1].to_i
            elsif (5..(5+n-1)).include?(ln)
                _, p, w, _ = line.split(',').map(&:to_i)
                weights << w
                profits << p
            elsif ln == 5+n
                return n, c, weights, profits
            end

            ln += 1
        end
    end
end

def print_instance(filename, n, c, w, p, pi)
    str = "#{n} #{c}\n"
    0.upto(n-1) do |i|
        str += "#{w[i]} #{p[i]} #{pi[i].round(3)}\n"
    end
    File.write(filename, str)
end

def generate_all_type1
    [100, 500, 1000, 5000].each do |n|
        [0.1, 0.2, 0.5].each do |b|
            1.upto(10) do |inst_n|
                w, p, pi, c = generate_type1(n, b)
                filename = "type1-#{n}-0-#{b}-#{inst_n}.txt"
                filename = File.join('generated-instances', filename)
                print_instance(filename, n, c, w, p, pi)
            end
        end
    end
end

def generate_all_type2
    [100, 500, 1000, 5000].each do |n|
        [0.1, 0.2, 0.5].each do |b|
            1.upto(10) do |inst_n|
                w, p, pi, c = generate_type2(n, b)
                filename = "type2-#{n}-0-#{b}-#{inst_n}.txt"
                filename = File.join('generated-instances', filename)
                print_instance(filename, n, c, w, p, pi)
            end
        end
    end
end

def generate_all_type3
    [100, 500, 1000, 5000].each do |n|
        [0.1, 0.2, 0.5].each do |b|
            1.upto(10) do |inst_n|
                w, p, pi, c = generate_type3(n, b)
                filename = "type3-#{n}-0-#{b}-#{inst_n}.txt"
                filename = File.join('generated-instances', filename)
                print_instance(filename, n, c, w, p, pi)
            end
        end
    end
end

def generate_all_type4
    [100, 500, 1000, 5000].each do |n|
        [0.1, 0.2, 0.5].each do |b|
            1.upto(10) do |inst_n|
                w, p, pi, c = generate_type4(n, b)
                filename = "type4-#{n}-0-#{b}-#{inst_n}.txt"
                filename = File.join('generated-instances', filename)
                print_instance(filename, n, c, w, p, pi)
            end
        end
    end
end

def generate_all_type5
    c = 10000

    [100, 500, 1000, 5000].each do |n|
        [1.6, 2.0, 2.4].each do |k|
            1.upto(10) do |inst_n|
                w, p, pi = generate_type5(n, k, c)
                filename = "type5-#{n}-#{k}-0-#{inst_n}.txt"
                filename = File.join('generated-instances', filename)
                print_instance(filename, n, c, w, p, pi)
            end
        end
    end
end

generate_all_type1
generate_all_type2
generate_all_type3
generate_all_type4
generate_all_type5