# adaptive-hermite-refinement

## Open Questions (assignee)
1. Do higher level moments have more spatial dimensions? (Alex)
2. Do we save intermediate moment values or only care about the final result? (Alex) Ans: If intermediate means at different timesteps, then yes, probably as a data file every n timesteps, controlled by user. (should not just be stored as giant array, that would be ridiculous)


## Governing Equations

**TODO: Vincent**:
- Hermite moment equation (7-12 in the paper)
- Spectral version
- Discretized version
  - just discretize the spatial dimension to compute the time derivative. It seems straightforward to take the time derivative and numerically integrate that with some scheme (Runge-Kutta or something else)

## Pseudocode

```
# setup grid
nMoments = Int[kx][ky] # how many moments at each point
moments, tmp_moments = Array{Moment}[kx][ky] # TODO moment dimensionality

# initialize, potentially via Fourier transform of input

# time-integration
for t in 1..T begin
  # don't refine every time ideally
  if t % REFINE_INTERVAL == 0 begin
    # for each point in space, check condition and update nMoments
    for kx in 1..KX begin
      for ky in 1..KY begin
        nMoments[kx][ky] = refineMoments(nMoments[ky][kx], moments[ky][kx])
      end
    end
  end
   
  for kx in 1..KX begin
    for ky in 1..KY begin
      # TODO do the first two moments
      for m in 3..nMoments[kx][ky]
        tmp_moments[kx][ky][m] = f(moments) # TODO ...
      end
    end
  end
  
  moments = tmp_moments
  # plot/output?
end
