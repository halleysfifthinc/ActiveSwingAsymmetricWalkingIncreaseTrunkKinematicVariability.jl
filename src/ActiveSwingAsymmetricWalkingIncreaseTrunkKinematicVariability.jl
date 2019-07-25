module ActiveSwingAsymmetricWalkingIncreaseTrunkKinematicVariability

using Biomechanics, MAT, Statistics, StaticArrays, DSP

export readsegment,
       readtrials,
       analyzetrial

export DSSteadyState,
       SteadyStateSeg

abstract type DSSteadyState <: AbstractDataSource end

struct SteadyStateSeg <: DSSteadyState
    events::Dict{Symbol,Vector{Float64}}
    data::Dict{Symbol,Matrix{Float64}}
end

SteadyStateSeg() = SteadyStateSeg(Dict{Symbol,Matrix}(), Dict{Symbol,AbstractVector}())

function readsegment(DSData::Type{<:AbstractDataSource},
                     trial::Trial{<:AbstractDataSource},
                     st::Float64 = 0.0,
                     en::Float64 = Inf;
                     events::Vector = [],
                     ts::Vector = [],
                     fs = 100)
    isfile(trial.path) || throw(ArgumentError("trial $(trial.path) exist"))
    st >= 0.0 || throw(ArgumentError("start time must be positive"))
    st <= en || throw(ArgumentError("end time must be greater than start time"))

    evnames = Symbol.(events)
    file = matopen(trial.path)

    revents = Dict{Symbol,Vector{Float64}}()
    for e in events
        if exists(file, string(e))
            syme = Symbol(e)
            tmp = read(file, string(e))[1]
            if tmp isa AbstractArray
                revents[syme] = vec(tmp)
            else
                revents[syme] = [tmp]
            end
            strt = findfirst(x -> x >= st, revents[syme])
            if strt == 0
                # @warn "no $e events during given start and end times"
                delete!(revents, syme)
                break
            end
            endi = findlast(x -> x <= en, revents[syme])
            revents[syme] = revents[syme][strt:endi] .- st .+ (1/fs) # Shift events to be index accurate for data subsection
        else
            # @warn "Requested event $e does not exist in source data"
        end
    end

    data = Dict{Symbol,Matrix{Float64}}()
    for t in ts
        if exists(file, string(t))
            symt = Symbol(t)
            data[symt] = read(file, string(t))[1]
            len = size(data[symt], 1)
            strti = round(Int, st*fs)
            endi = en == Inf ? len : min(len, round(Int, en*fs))
            data[symt] = data[symt][strti:endi, :]
        else
            @warn "Requested time series $t does not exist in source data"
        end
    end

    return Segment(trial, Dict{Symbol,Any}(), DSData(revents, data))
end

const fs = 100

function readtrials(datadir::String)
    sstrials = Vector{Trial}()

    fp = abspath(datadir)

    trials = readdir(datadir)
    filter!(trial -> endswith(trial, "mat"), trials)

    for trial in trials
        conds = Dict{Symbol, Symbol}()
        m = match(r"S(?<sub>\d{1,2})_(?<arms>(held|norm|active))_(?<sym>(sym|asym))", trial)
        sub = tryparse(Int, m[:sub])
        conds[:sym] = Symbol(m[:sym])
        conds[:arms] = Symbol(m[:arms])

        push!(sstrials, Trial{DSSteadyState}(sub, splitext(basename(trial))[1], joinpath(fp, trial), conds))
    end

    return sstrials
end

function analyzetrial(trial, numstrides)
    cols = [ "TrunkLinVel",
             "TrunkAngVel",
             "RFootPos",
             "LFootPos",
             "ModelAngMmntm" ]
    seg = readsegment(SteadyStateSeg, trial, 25.0; events=[:RFC, :LFC], ts=cols)
    insuffstrerr = "insufficient number of strides in trial $(trial.name), subject $(trial.subject)"
    length(seg.data.events[:RFC]) < numstrides+1 && throw(ArgumentError(insuffstrerr))

    results = Dict{Symbol,Any}()

    ########################################
    ## Interleaving steps
    ########################################

    # Gather rfc and lfc info into a NamedTuple, with named fields `time`, `pos`, and `leg`, (`:RFC` or `:LFC`)
    rfc = [ (time=seg.data.events[:RFC][i],
             pos=SVector{2}(seg.data.data[:RFootPos][
                            round(Int, seg.data.events[:RFC][i]*fs)
                            ,1:2]),
             leg=:RFC) for i in 1:(numstrides+1) ]
    if seg.data.events[:LFC][1] < seg.data.events[:RFC][1]
        lfc = [ (time=seg.data.events[:LFC][i],
                 pos=SVector{2}(seg.data.data[:LFootPos][
                                round(Int, seg.data.events[:LFC][i]*fs)
                                ,1:2]),
                 leg=:LFC) for i in 2:(numstrides+2) ]
    else
        lfc = [ (time=seg.data.events[:LFC][i],
                 pos=SVector{2}(seg.data.data[:LFootPos][
                                round(Int, seg.data.events[:LFC][i]*fs)
                                ,1:2]),
                 leg=:LFC) for i in 1:(numstrides+1) ]
    end

    # Sort all steps and check that they are alternating--no double steps--and that they are
    # in the expected locations:
    # Odd elements need to be RFC for later assumptions made on whether the odd/even elements
    # of the diff are right or left steps
    steps = sort([rfc; lfc], by=(x -> x.time))

    doublesteps = 0
    while any(x -> x.leg === :LFC, steps[isodd.(axes(steps, 1))])
        badstep = min(something(findfirst(x -> x.leg === :LFC, steps[1:2:end]), typemax(Int)÷2)*2-1,
                      something(findfirst(x -> x.leg === :RFC, steps[2:2:end]), typemax(Int)÷2)*2)
        steps = deleteat!(steps, badstep)
        doublesteps += 1
    end
    all(x -> x.leg === :RFC, steps[isodd.(axes(steps, 1))]) || throw(DomainError("Odd elements aren't all RFC's"))
    all(x -> x.leg === :LFC, steps[iseven.(axes(steps, 1))]) || throw(DomainError("Even elements aren't all LFC's"))
    if doublesteps != 0
        @warn "Double stepping occured" trial, doublesteps
    end

    ########################################
    ## Gait variability
    ########################################

    # After `diff`ing, odd elements will be left steps, even elements will be right steps
    steptimes = diff(getindex.(steps, :time))

    if trial.conds[:sym] == :asym
        # The speed of the treadmill belt for the planted foot is added to the distance 
        # between successive steps e.g. For the left step, the right foot is planted,
        # therefore the effective AP distance traveled by the left foot is speed of the right
        # treadmill belt plus the absolute distance between the right and left footstrikes
        stepcoords = diff(getindex.(steps, :pos)) .+
                        ( isodd(i) ? SVector{2}(0.0, 0.96*steptimes[i]) : SVector{2}(0.0, 1.2*steptimes[i])
                         for i in 1:(length(steps)-1) )
    else
        stepcoords = diff(getindex.(steps, :pos)) .+
                        ( SVector{2}(0.0, 1.2*steptimes[i]) for i in 1:(length(steps)-1) )
    end

    # Spatial and temporal step descriptives
    left_stepcoords = stepcoords[isodd.(axes(stepcoords,1))]
    right_stepcoords = stepcoords[iseven.(axes(stepcoords,1))]
    len = min(length(left_stepcoords), length(right_stepcoords))
    
    left_stepcoords = left_stepcoords[1:len]
    right_stepcoords = right_stepcoords[1:len]
    
    results[:left_steplength] = mean(abs.(getindex.(left_stepcoords,2)))
    results[:SD_left_steplength] = std(abs.(getindex.(left_stepcoords,2)))

    results[:right_steplength] = mean(abs.(getindex.(right_stepcoords,2)))
    results[:SD_right_steplength] = std(abs.(getindex.(right_stepcoords,2)))

    results[:left_steptime] = mean(steptimes[isodd.(axes(steptimes,1))])
    results[:right_steptime] = mean(steptimes[iseven.(axes(steptimes,1))])
    results[:SD_left_steptime] = std(steptimes[isodd.(axes(steptimes,1))])
    results[:SD_right_steptime] = std(steptimes[iseven.(axes(steptimes,1))])

    # Step width isn't confounded by induced asymmetric gait
    results[:stepwidth] = mean(abs.(getindex.(stepcoords, 1)))
    results[:SD_stepwidth] = std(abs.(getindex.(stepcoords, 1)))

    ########################################
    ## LV and AV reductions
    ########################################

    LV = seg.data.data[:TrunkLinVel][1:end-1,:]
    AV = seg.data.data[:TrunkAngVel][1:end-1,:]
    
    stridergs = [ round(Int, seg.data.events[:RFC][i]*fs):(round(Int, seg.data.events[:RFC][i+1]*fs)-1) 
                    for i in 1:numstrides ]
    
    results[:lvmean] = mean([ vec(mean(view(LV, rg, :); dims=1)) for rg in stridergs ])
    results[:avmean] = mean([ vec(mean(view(AV, rg, :); dims=1)) for rg in stridergs ])
    
    results[:lvstd] = mean([ vec(std(view(LV, rg, :); dims=1)) for rg in stridergs ])
    results[:avstd] = mean([ vec(std(view(AV, rg, :); dims=1)) for rg in stridergs ])
    
    results[:lvmax] = mean([ vec(maximum(view(LV, rg, :); dims=1)) for rg in stridergs ])
    results[:avmax] = mean([ vec(maximum(view(AV, rg, :); dims=1)) for rg in stridergs ])
        
    ########################################
    ## WBAM reductions
    ########################################

    AM = seg.data.data[:ModelAngMmntm][3:end-2,:]
    
    stridergs .= [ round(Int, seg.data.events[:RFC][i]*fs):(round(Int, seg.data.events[:RFC][i+1]*fs)-1) .- 2 # Adjust for the 2 removed samples
                    for i in 1:numstrides ]

    results[:WBAM_mean] = mean([ vec(mean(view(AM, rg, :); dims=1)) for rg in stridergs ])
    results[:WBAM_std] = mean([ vec(std(view(AM, rg, :); dims=1)) for rg in stridergs ])

    return AnalyzedSegment(seg, results)
end

end # module
