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

generate_all_uv