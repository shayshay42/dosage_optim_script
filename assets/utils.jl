using Plots, LinearAlgebra, Statistics

"""
    trapezoidal_rule(t, y)

Calculate the area under a curve defined by 'x' and 'y' coordinates using the trapezoidal rule. This numerical integration method works by approximating the region under the curve as a series of trapezoids and then calculating the sum of their areas.

# Arguments
- 't::Vector': A vector of 'x' values (time points) in ascending order.
- 'y::Vector': A vector of 'y' values (concentration levels) corresponding to each 'x' value.

# Returns
- 'Float': The approximate area under the curve, as calculated by the trapezoidal rule.

# Example
'''julia
julia> t = 0:0.1:1;
julia> y = t.^2;
julia> trapezoidal_rule(t, y)
0.33000000000000007
"""
function trapezoidal_rule(t, y)
  return 0.5*sum((t[2:end] - t[1:end-1]) .* (y[1:end-1] + y[2:end]))
end

"""
"""
function relu(x)
  return max(0,x)
end

"""
    spaced_list(p, n, m, b=1)

Create a list of integers where 'n' numbers are sequentially appended 
followed by a jump of 'm' numbers. This pattern is repeated until 
the end number 'p' is reached or surpassed. 

Optionally, the start of the sequence can be adjusted from the 
default of '1' with the 'b' parameter.

# Arguments
- 'p::Integer': The final number of the sequence. The function will stop 
  adding numbers to the list once this number is reached or surpassed.
- 'n::Integer': The number of sequential integers to append to the list 
  at a time.
- 'm::Integer': The number of integers to skip in the sequence after 
  each set of 'n' integers is added.
- 'b::Integer': (optional) The beginning number of the sequence. Default 
  is '1'.

# Returns
- 'Array{Integer}': An array of integers that follows the specified 
  sequential pattern.

# Example
'''julia
julia> spaced_list(20, 2, 3)
[1, 2, 6, 7, 11, 12, 16, 17]
"""
function spaced_list(p, n, m, b=1)
  # Can be optimized...
  
  spaced_list = []
  #iniialize counter
  counter = b
  #range of integers
  while counter <= p
    #append n integers spaced by 1 to the list
    for i in 1:n
      spaced_list = [spaced_list; counter]
      counter += 1

      if counter > p
        break
      end
    end
    #add m to create the jump
    counter += m
  end
  return spaced_list
end

"""
    create_callback()

Create and return a callback function for an optimization routine and a list to store loss values.

The returned callback function prints the epoch number and loss value at each step of the optimization,
and also stores each loss value in the returned list. The callback is intended to be used with an 
optimization routine from the Optimization.jl package.

# Returns
- 'track_progress::Function': the callback function that takes 'iter' (the current epoch number) and 
  'x' (the current parameter values) as arguments, computes the loss, and stores it in 'loss_values'.
- 'loss_values::Vector': a list to hold the loss values computed during the optimization.

# Example
'''julia
callback, loss_values = create_callback()
options = Optimization.Options(callback = callback, maxiters = 50)
res = Optimization.solve(optprob, opt, options)
@show res.u
@show loss_values

"""
function create_callback()
  loss_values = []
  function track_progress(iter, x)
      current_loss = loss(x)

      println("Epoch: $iter, Loss: $current_loss")

      #save loss in local list
      push!(loss_values, current_loss)
      
      #continue the optimization
      return false  
  end
  return track_progress, loss_values
end