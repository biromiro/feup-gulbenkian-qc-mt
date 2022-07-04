// Programming Quantum Computers
//   by Eric Johnston, Nic Harrigan and Mercedes Gimeno-Segovia
//   O'Reilly Media

// To run this online, go to http://oreilly-qc.github.io?p=12-1
// Note: This sample may vary slightly from the text in the book,
// due to revisions or aesthetic tweaks.

// Special note: This implementation of Shor's algorithm is for
// illustration purposes, to help develop an intuition regarding
// what the algorithm does. It is not intended to be an optimal
// implementation on any specific QPU or simulation.

// Here are some values of N to try:
// 15, 21, 35, 39, 51, 55, 69, 77, 85, 87, 91, 93, 95, 111, 115, 117,
// 119, 123, 133, 155, 187, 203, 221, 247, 259, 287, 341, 451

// Larger numbers require more bits of precision.
// N = 15    precision_bits >= 4
// N = 21    precision_bits >= 5
// N = 35    precision_bits >= 6
// N = 123   precision_bits >= 7
// N = 341   precision_bits >= 8  time: about 6 seconds
// N = 451   precision_bits >= 9  time: about 23 seconds

var increase_speed = false; // switch drawing off to increase sim speed
var primes_list = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37,
               41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97];

var prime_square = new Map();
var prime_cube = new Map();

for (let prime in primes_list) {
    prime_square.set(prime*prime, prime);
    prime_cube.set(prime*prime*prime, prime);
}

    
function precision_bits(n) {
    if (n <= 15) return 4;
    if (n <= 21) return 5;
    if (n <= 35) return 6;
    if (n <= 123) return 7;
    if (n <= 341) return 8;
    if (n <= 451) return 9;
    return 10;
}

function execute_shor(val, a) {
    
    let result = null;
    let attempts = 20;

    while (result == null && attempts > 0) {
        result = Shor(val, precision_bits(val), a);
        attempts -= 1;
    }
    
    return result;
}

function factorize(result) {
    let factorized = false;
    
    while (!factorized) {
        factorized = true;
        let new_result = [];
        
        for (let coprime of result) {
    		let prime = parseInt(prime_square.get(coprime));
    		if (prime) {
    		      factorized = false;
                  new_result = new_result.concat([prime, prime]);
    		      continue;
    		}
    		prime = parseInt(prime_cube.get(coprime));
    		if (prime) {
    		      factorized = false;
                  new_result = new_result.concat([prime, prime, prime]);
    		      continue;
    		}
            if (!primes_list.includes(coprime)) {
            	factorized = false;
            	let new_factors = null;
            	for (let a of [2] + primes_list) {
                    if (new_factors != null || a > coprime) break;
            	    new_factors = execute_shor(coprime, a);
            	}
            	if (new_factors != null) {
            		new_result = new_result.concat(new_factors);
            	} else new_result.push(coprime);
            } else new_result.push(coprime);
        }
        
        result = new_result;
        result.sort();
    }
    
    return result;
}

function shor_sample()
{

    for (let prime of primes_list) {
        for (var q = 15; q <= 35; q += 2) {
            let q2 = prime * q;

            if (primes_list.includes(q) || q < prime || q2 > 450) continue;
            
            let result1 = null, result2 = null, t1, t2, t3, t4;
            
            for (let a of [2] + primes_list) {
                if (result1 != null || a > q) break;
                t1 = performance.now();
                result1 = execute_shor(q, a);
                t2 =  performance.now();
            }
            
            result1 = factorize(result1);
           
            for (let a of [2] + primes_list) {
                if (result2 != null || a > q2) break;
                t3 = performance.now();
                result2 = execute_shor(q2, a);
                t4 =  performance.now();
            }
            
            result2 = factorize(result2);
            
            const q1_time = t2 - t1;
            const q2_time = t4 - t3;

            if (result1 != null && result2 != null)
                qc.print('With prime ' + prime + ':\nFactors of q(' + q + ')(' + q1_time +'ms): ' + result1 + '\n Factors of q2(' + q2 + ')(' + q2_time +'ms): ' + result2 + '\n');
            else if (result1 != null)
                qc.print('With prime ' + prime + ':\nFactors of q(' + q + ')(' + q1_time +'ms): ' + result1 + '\n Couldn\'t factor q2(' + q2 + ')\n');
            else if (result2 != null)
                qc.print('With prime ' + prime + ':\nCouldn\'t factor q(' + q + ')\n Factors of q2(' + q2 + ')(' + q2_time +'ms): ' + result2 + '\n');
            else
                qc.print('With prime ' + prime + ':\nFailure: No non-trivial factors were found (' + q + ',' + q2 + ')\n')
            
        }
    }

}

function Shor(N, precision_bits, coprime)
{
    var repeat_period = ShorQPU(N, precision_bits, coprime); // quantum part
    var factors = ShorLogic(N, repeat_period, coprime);      // classical part
    return check_result(N, factors);
}

function gcd(a, b)
{
    // return the greatest common divisor of a,b
    while (b) {
      var m = a % b;
      a = b;
      b = m;
    }
    return a;
}

function check_result(N, factor_candidates)
{
    for (var i = 0; i < factor_candidates.length; ++i)
    {
        var factors = factor_candidates[i];
        if (factors[0] * factors[1] == N)
        {
            if (factors[0] != 1 && factors[1] != 1)
            {
                // Success!
                return factors;
            }
        }
    }
    // Failure
    return null;
}

function ShorLogic(N, repeat_period_candidates, coprime)
{
    //qc.print('Repeat period candidates: '+repeat_period_candidates+'\n');
    factor_candidates = [];
    for (var i = 0; i < repeat_period_candidates.length; ++i)
    {
        var repeat_period = repeat_period_candidates[i];
    // Given the repeat period, find the actual factors
        var ar2 = Math.pow(coprime, repeat_period / 2.0);
        var factor1 = gcd(N, ar2 - 1);
        var factor2 = gcd(N, ar2 + 1);
        factor_candidates.push([factor1, factor2]);
    }
    return factor_candidates;
}

function setup_speed()
{
    if (increase_speed)
    {
        qc.disableRecording();
        qc.disableAnimation();
    }
    else
    {
        qc.enableRecording();
        qc.enableAnimation();
    }
}

function ShorQPU(N, precision_bits, coprime)
{
    setup_speed();
    // Quantum part of Shor's algorithm
    // For this implementation, the coprime must be 2.
    coprime = 2;

    // For some numbers (like 15 and 21) the "mod" in a^xmod(N)
    // is not needed, because a^x wraps neatly around. This makes the
    // code simpler, and much easier to follow.
    if (N == 15 || N == 21)
        return ShorQPU_WithoutModulo(N, precision_bits, coprime)
    else
        return ShorQPU_WithModulo(N, precision_bits, coprime)
}

// In case our QPU read returns a "signed" negative value,
// convert it to unsigned.
function read_unsigned(qreg)
{
    var value = qreg.read();
    return value & ((1 << qreg.numBits) - 1);
}

// This is the short/simple version of ShorQPU() where we can perform a^x and
// don't need to be concerned with performing a quantum int modulus.
function ShorQPU_WithoutModulo(N, precision_bits, coprime)
{
    var N_bits = 1;
    while ((1 << N_bits) < N)
    N_bits++;
    if (N != 15) // For this implementation, numbers other than 15 need an extra bit
    N_bits++;
    var total_bits = N_bits + precision_bits;

    // Set up the QPU and the working registers
    qc.reset(total_bits);
    var num = qint.new(N_bits, 'work');
    var precision = qint.new(precision_bits, 'precision');

    qc.label('init');
    num.write(1);
    precision.write(0);
    precision.hadamard();

    // Perform 2^x for all possible values of x in superposition
    for (var iter = 0; iter < precision_bits; ++iter)
    {
        qc.label('iter ' + iter);
        var num_shifts = 1 << iter;
        var condition = precision.bits(num_shifts);
        num_shifts %= num.numBits;
        num.rollLeft(num_shifts, condition);
    }
    // Perform the QFT
    qc.label('QFT');
    precision.QFT();
    qc.label('');

    var read_result = read_unsigned(precision);
    //qc.print('QPU read result: '+read_result+'\n')
    var repeat_period_candidates = estimate_num_spikes(read_result, 1 << precision_bits);

    return repeat_period_candidates;
}

// This is the complicated version of ShorQPU() where we DO
// need to be concerned with performing a quantum int modulus.
// That's a complicated operation, and it also requires us to
// do the shifts one at a time.
function ShorQPU_WithModulo(N, precision_bits, coprime)
{
    var scratch = null;
    var max_value = 1;
    var mod_engaged = false;

    var N_bits = 1;
    var scratch_bits = 0;
    while ((1 << N_bits) < N)
        N_bits++;
    if (N != 15) // For this implementation, numbers other than 15 need an extra bit
        N_bits++;
    scratch_bits = 1;
    var total_bits = N_bits + precision_bits + scratch_bits;

    // Set up the QPU and the working registers
    qc.reset(total_bits);
    var num = qint.new(N_bits, 'work');
    var precision = qint.new(precision_bits, 'precision');
    var scratch = qint.new(1, 'scratch');

    qc.label('init');
    num.write(1);
    precision.write(0);
    precision.hadamard();
    scratch.write(0);

    var N_sign_bit_place = 1 << (N_bits - 1);
    var N_sign_bit = num.bits(N_sign_bit_place);
    for (var iter = 0; iter < precision_bits; ++iter)
    {
        var condition = precision.bits(1 << iter);
        var N_sign_bit_with_condition = num.bits(N_sign_bit_place);
        N_sign_bit_with_condition.orEquals(condition);

        var shifts = 1 << iter;
        for (var sh = 0; sh < shifts; ++sh)
        {
            qc.label('num *= coprime');
            num.rollLeft(1, condition);   // Multiply by the coprime
            max_value <<= 1;
            if (max_value >= N)
                mod_engaged = true;
            if (mod_engaged)
            {
                qc.label('modulo N');
                var wrap_mask = scratch.bits();
                var wrap_mask_with_condition = scratch.bits();
                wrap_mask_with_condition.orEquals(condition);

                // Here's the modulo code.
                num.subtract(N, condition); // subtract N, causing this to go negative if we HAVEN'T wrapped.
                scratch.cnot(N_sign_bit_with_condition); // Skim off the sign bit
                num.add(N, wrap_mask_with_condition); // If we went negative, undo the subtraction.
                num.not(1);
                scratch.cnot(num, 1, condition); // If it's odd, then we wrapped, so clear the wrap bit
                num.not(1);
            }
        }
    }
    qc.label('QFT');
    precision.QFT();
    qc.label('');

    var read_result = read_unsigned(precision);
    //qc.print('QPU read result: '+read_result+'\n')
    var repeat_period_candidates = estimate_num_spikes(read_result, 1 << precision_bits);

    return repeat_period_candidates;
}

function estimate_num_spikes(spike, range)
{
    if (spike < range / 2)
        spike = range - spike;
    var best_error = 1.0;
    var e0 = 0;
    var e1 = 0;
    var e2 = 0;
    var actual = spike / range;
    var candidates = []
    for (var denom = 1.0; denom < spike; ++denom)
    {
        var numerator = Math.round(denom * actual);
        var estimated = numerator / denom;
        var error = Math.abs(estimated - actual);
        e0 = e1;
        e1 = e2;
        e2 = error;
        // Look for a local minimum which beats our
        // current best error
        if (e1 <= best_error && e1 < e0 && e1 < e2)
        {
            var repeat_period = denom - 1;
            candidates.push(repeat_period);
            best_error = e1;
        }
    }
    return candidates;
}

shor_sample();
