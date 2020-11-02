require 'simple-random'

def generate_uu(n, k, c, b)
    w = (k * c / n).floor.to_i
    weights = Array.new(n) {w}

    r = SimpleRandom.new
    profits = Array.new(n) {[1, (w + r.normal(0, w/8.0)).ceil].max}.map(&:to_i)

    pprod = 0.5 * profits.min
    probabilities = Array.new(n) do |i|
        (rand() < b) ? pprod / profits[i] : 1
    end.map{|pi| pi.round(6)}

    return weights, profits, probabilities
end

def generate_uv(n, k, c, b)
    w = (k * c / n).floor.to_i
    weights = Array.new(n) {w}

    r = SimpleRandom.new
    profits = Array.new(n) {[1, (w + r.normal(0, w/8.0)).ceil].max}.map(&:to_i)

    probabilities = Array.new(n) do
        (rand() < b) ? (1.0 - r.beta(2, 10)) : 1
    end.map{|pi| pi.round(6)}

    return weights, profits, probabilities
end

def generate_vu(n, b)
    f = Dir.glob(File.join('from-pisinger-hard', "knapPI_*_#{n}_*.csv")).to_a.sample
    i = rand(1..100)
    nn, capacity, weights, profits = get_pis_instance(f, i)

    raise "Wrong n from Pisinger instance (exp #{n}, got #{nn})!" if n != nn

    pprod = 0.5 * profits.min
    probabilities = Array.new(n) do |i|
        (rand() < b) ? pprod / profits[i] : 1
    end.map{|pi| pi.round(6)}

    return weights, profits, probabilities, capacity
end

def generate_vv(n, b)
    f = Dir.glob(File.join('from-pisinger-hard', "knapPI_*_#{n}_*.csv")).to_a.sample
    i = rand(1..100)
    nn, capacity, weights, profits = get_pis_instance(f, i)

    raise "Wrong n from Pisinger instance (exp #{n}, got #{nn})!" if n != nn

    r = SimpleRandom.new
    probabilities = Array.new(n) do
        (rand() < b) ? (1.0 - r.beta(2, 10)) : 1
    end.map{|pi| pi.round(6)}

    return weights, profits, probabilities, capacity
end

def generate_c(n, b)
    f = Dir.glob(File.join('from-pisinger-hard', "knapPI_*_#{n}_*.csv")).to_a.sample
    i = rand(1..100)
    nn, capacity, weights, profits = get_pis_instance(f, i)

    raise "Wrong n from Pisinger instance (exp #{n}, got #{nn})!" if n != nn

    r = SimpleRandom.new
    probabilities = Array.new(n) do
        (rand() < b) ? (1.0 - r.beta(2, 10)) : 1
    end.map{|pi| pi.round(6)}

    profits = profits.each_with_index.map do |p, id|
        p * (1 + 2 * (1 - probabilities[id]))
    end

    return weights, profits, probabilities, capacity
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
        str += "#{w[i]} #{p[i]} #{pi[i]}\n"
    end
    File.write(filename, str)
end

def generate_all_uu
    c = 100000

    [100, 500, 1000, 5000].each do |n|
        [1.6, 1.8, 2.0, 2.2, 2.4].each do |k|
            bs = [0.05, 0.1, 0.2, 0.5]
            bs << 0.01 if n >= 1000

            bs.each do |b|
                1.upto(5) do |inst_n|
                    w, p, pi = generate_uu(n, k, c, b)
                    filename = "uu-#{n}-#{k}-#{b}-#{inst_n}.txt"
                    filename = File.join('generated-instances', filename)
                    print_instance(filename, n, c, w, p, pi)
                end
            end
        end
    end
end

def generate_all_uv
    c = 100000

    [100, 500, 1000, 5000].each do |n|
        [1.6, 1.8, 2.0, 2.2, 2.4].each do |k|
            bs = [0.05, 0.1, 0.2, 0.5]
            bs << 0.01 if n >= 1000

            bs.each do |b|
                1.upto(5) do |inst_n|
                    w, p, pi = generate_uv(n, k, c, b)
                    filename = "uv-#{n}-#{k}-#{b}-#{inst_n}.txt"
                    filename = File.join('generated-instances', filename)
                    print_instance(filename, n, c, w, p, pi)
                end
            end
        end
    end
end

def generate_all_vu
    [100, 500, 1000, 5000].each do |n|
        bs = [0.05, 0.1, 0.2, 0.5]
        bs << 0.01 if n >= 1000

        bs.each do |b|
            1.upto(25) do |inst_n|
                w, p, pi, c = generate_vu(n, b)
                filename = "vu-#{n}-0-#{b}-#{inst_n}.txt"
                filename = File.join('generated-instances', filename)
                print_instance(filename, n, c, w, p, pi)
            end
        end
    end
end

def generate_all_vv
    [100, 500, 1000, 5000].each do |n|
        bs = [0.05, 0.1, 0.2, 0.5]
        bs << 0.01 if n >= 1000

        bs.each do |b|
            1.upto(25) do |inst_n|
                w, p, pi, c = generate_vv(n, b)
                filename = "vv-#{n}-0-#{b}-#{inst_n}.txt"
                filename = File.join('generated-instances', filename)
                print_instance(filename, n, c, w, p, pi)
            end
        end
    end
end

def generate_all_c
    [100, 500, 1000, 5000].each do |n|
        bs = [0.05, 0.1, 0.2, 0.5]
        bs << 0.01 if n >= 1000

        bs.each do |b|
            1.upto(25) do |inst_n|
                w, p, pi, c = generate_vv(n, b)
                filename = "c-#{n}-0-#{b}-#{inst_n}.txt"
                filename = File.join('generated-instances', filename)
                print_instance(filename, n, c, w, p, pi)
            end
        end
    end
end
