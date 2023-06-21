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
  iter = 1
  anim = Animation()
  function track_progress(x, current_loss)
      println("Epoch: $iter, Loss: $current_loss")
      iter += 1
      #save loss in local list
      push!(loss_values, current_loss)
      frame(anim, plotter(sim(x)))
      #continue the optimization
      return false  
  end
  return track_progress, loss_values, iter, anim
end

"""
    plotter(sol)

This function generates two plots with time on the x-axis, and the values of doses and cells on the y-axes, respectively.
The data is dynamically adjusted according to the ranges in the provided data `sol`.

# Arguments
- `sol::Array`: A 2D array where the first row corresponds to the cell values, the fifth and eighth rows correspond to the dose values. The time values are stored in `sol.t`.

# Returns
- `Plots.Plot`: A plot object generated by the Plots package, showing two subplots stacked vertically. 

The first subplot ("Zoom for CSF concentrations") has the dose amounts for TMZ in CSF and GDC in the periphery plotted on the left y-axis, and the cell volume plotted on the right y-axis. A line of annotation indicating the final volume is also included.

The second subplot ("Dosing Regimen") has the dose amounts for AbsTMZ, PlaTMZ, CSFTMZ, AbsGDC, PlaGDC, and PeriphGDC plotted on the left y-axis, and the cell volume plotted on the right y-axis.

# Example
```julia
julia> plotter(sol)
"""
function plotter(sol)
  c = sol[1,:]
  fv = round(c[end], digits=4)
  cells = Array(sol)[1:1,:]
  doses = Array(sol)[[5,7],:]
  time = sol.t

  dose_min = minimum(transpose(doses))
  dose_max = maximum(transpose(doses))

  cell_min = minimum(transpose(cells))
  cell_max = maximum(transpose(cells))

  p1 = plot(time./24, transpose(doses), ylims=(dose_min * 0.9, dose_max * 1.1), color = ["#9bcfa0" "#d6818f"], title="Dosing Simulation", linewidth=1, xaxis="time (days)", label=["TMZ in CSF" "Drug in Plasma"], ylabel="Dose Amount (mg)", legend=:topleft)
  p2 = twinx(p1)
  plot!(p2, time./24, transpose(cells), ylims=(cell_min * 0.9, cell_max * 1.1), color = "#81b1d6", linewidth=3, linestyle = :dash, label="C", ylabel="Tumor Volume (ml)", legend=:topright)

  mid_x = (maximum(time) - minimum(time))/2 + minimum(time)
  mid_y = (maximum(cells) - minimum(cells))/2 + minimum(cells)

  annotate!(mid_x, mid_y, text("\n Final Volume: $fv", :black, :left, 10))
  return p1
end

function sim(x)
  q = [ode_params;x]
  tmp_prob = remake(prob, p=q)
  tmp_sol = solve(tmp_prob, p=q, callback=hit)
  return tmp_sol
end