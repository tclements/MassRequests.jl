using SeisRequests
using CSV
using DataFrames
using Dates
import DataStructures: OrderedDict

const MBool = Union{Missing,Bool}
const MDateTime = Union{Missing,DateTime}
const MFloat = Union{Missing,Float64}
const MInt = Union{Missing,Int}
const MString = Union{Missing,String}

"Available processing filters and allowable types for them"
const IRISTimeSeries_PROCESSING_FIELDS = Dict(
    :taper => Real,
    :taper_type => AbstractString,
    :envelope => Bool,
    :lpfilter => Real,
    :hpfilter => Real,
    :bpfilter => Tuple{Number,Number},
    :demean => Bool,
    :scale => Union{Real,AbstractString},
    :divscale => Real,
    :correct => Bool,
    :freqlimits => Tuple{Number,Number,Number,Number},
    :autolimits => Union{AbstractString,Tuple{Number,Number}},
    :units => AbstractString,
    :diff => Bool,
    :int => Bool,
    :decimate => Real)

const IRISTimeSeries_OUTPUT_TYPES = ("ascii", "ascii2", "ascii", "geocsv", "geocsv.tspair",
"geocsv.slist", "audio", "miniseed", "plot", "saca", "sacbb", "sacbl")

struct IRISMassRequest
  # Obspy Domain & Restriction class
  network::String
  station::String
  location::String
  channel::String
  starttime::DateTime
  endtime::DateTime
  duration::MFloat
  output::String
  process::OrderedDict{Symbol,Any}
  minlatitude::MFloat
  maxlatitude::MFloat
  minlongitude::MFloat
  maxlongitude::MFloat
  longitude::MFloat
  latitude::MFloat
  minradius::MFloat
  maxradius::MFloat
  startbefore::MDateTime
  startafter::MDateTime
  endbefore::MDateTime
  endafter::MDateTime
  reject_channels_with_gaps::MBool
  minimum_length::MFloat
  minimum_interstation_distance_in_m::MInt


"""


"""
  function IRISMassRequest(network, station, location, channel, starttime, endtime, duration,
                           output, process, minlatitude, maxlatitude, minlongitude, maxlongitude,
                            longitude, latitude, minradius, maxradius, startbefore, startafter,
                            endbefore, endafter, reject_channels_with_gaps, minimum_length,
                            minimum_interstation_distance_in_m)
        location == "  " && (location = "--")
        !ismissing(starttime) && !ismissing(endtime) && starttime > endtime &&
            throw(ArgumentError("`starttime` must be before endtime"))
        coalesce(duration, 1) > 0 || throw(ArgumentError("`duration` must be positive"))
        0 <= get(process, :taper, 0.1) <= 1 || throw(ArgumentError("`taper` length must be between 0 and 1"))
        # Check types for processing
        for (f, T) in IRISTimeSeries_PROCESSING_FIELDS
            if haskey(process, f)
                try
                    process[f] = convert(T, process[f])
                catch err
                    err isa MethodError || err isa InexactError &&
                        throw(ArgumentError("field `$f` must be of type $T"))
                    rethrow(err)
                end
            end
        end
        if haskey(process, :taper_type)
            haskey(process, :taper) || throw(ArgumentError("`taper_type` must be specficied with taper"))
            taper_type = uppercase(process[:taper_type])
            taper_type in ("HANNING", "HAMMING", "COSINE") ||
            throw(ArgumentError("taper_type must be one of \"HANNING\", \"HAMMING\" or \"COSINE\""))
        end
        haskey(process, :scale) && haskey(process, :divscale) &&
            throw(ArgumentError("Cannot specify both `scale` and `divscale`"))
        if haskey(process, :scale)
            if process[:scale] isa AbstractString
                process[:scale] = uppercase(process[:scale])
                process[:scale] == "AUTO" || throw(ArgumentError("`scale` must be a float or \"AUTO\""))
            end
        end
        if haskey(process, :units)
            !haskey(process, :correct) || !process[:correct] &&
                throw(ArgumentError("`units` can only be used when correct is true"))
            process[:units] in ("DIS", "VEL", "ACC", "DEF") ||
                throw(ArgumentError("`units` must be one of \"DIS\", \"VEL\", \"ACC\" or \"DEF\""))
        end
        haskey(process, :autolimits) && haskey(process, :freqlimits) &&
            throw(ArgumentError("Cannot specify both `autolimits` and `freqlimits`"))

        # Disallowed values
        output in IRISTimeSeries_OUTPUT_TYPES ||
            throw(ArgumentError("`output` must be one of $(IRISTimeSeries_OUTPUT_TYPES)"))

        # MassRequest values


        new(network, station, location, channel, starttime, endtime, duration,
            output, process, minlatitude, maxlatitude, minlongitude, maxlongitude,
            longitude, latitude, minradius, maxradius, startbefore, startafter,
            endbefore, endafter, reject_channels_with_gaps, minimum_length,
            minimum_interstation_distance_in_m)
    end
end

function IRISMassRequest(;
                        network=missing, station=missing, location=missing, channel=missing,
                        starttime=missing, endtime=missing, duration=missing,
                        output=missing, minlatitude=missing, maxlatitude=missing,
                        minlongitude=missing, maxlongitude=missing, longitude=missing,
                        latitude=missing, minradius=missing, maxradius=missing,
                        startbefore=missing, startafter=missing, endbefore=missing,
                        endafter=missing, reject_channels_with_gaps=missing,
                        minimum_length=missing, minimum_interstation_distance_in_m=missing,
                        kwargs...)
    any(ismissing, (network, station, location, channel, starttime, output)) &&
        throw(ArgumentError("network, station, location, channel, starttime " *
                            "and output must all be specified"))
    all(ismissing, (endtime, duration)) &&
        throw(ArgumentError("One of `endtime` or `duration` must be specified"))
    process = OrderedDict{Symbol,Any}()
    for (k, v) in kwargs
        k in keys(IRISTimeSeries_PROCESSING_FIELDS) ||
            throw(ArgumentError("Field `$k` is not a valid IRISTimeSeries processing field"))
        process[k] = v
    end
    IRISMassRequest(network, station, location, channel, starttime, endtime, duration,
                    lowercase(output), process, minlatitude, maxlatitude, minlongitude,
                    maxlongitude, longitude, latitude, minradius, maxradius, startbefore,
                    startafter, endbefore, endafter, reject_channels_with_gaps, minimum_length,
                   minimum_interstation_distance_in_m)
end

"""
  station_request(station,network,channel)


"""
function station_availability(mass_request::IRISMassRequest;
                         level::String="channel", format::String="text")

  # if ismissing(starttime)
  #   starttime = DateTime(1900,1,1)
  # end
  # if ismissing(endtime)
  #   endtime = DateTime(2599,12,31,23,59,59)
  # end

  request = FDSNStation(station = mass_request.station, network = mass_request.network,
                        channel = mass_request.channel,
                        level = level, format = format,
                        starttime = mass_request.starttime,
                        endtime = mass_request.endtime,
                        startbefore = mass_request.startbefore,
                        startafter = mass_request.startafter,
                        endbefore = mass_request.endbefore,
                        endafter = mass_request.endafter,
                        location = mass_request.location,
                        minlatitude = mass_request.minlatitude,
                        maxlatitude = mass_request.maxlatitude,
                        minlongitude = mass_request.minlongitude,
                        maxlongitude = mass_request.maxlongitude,
                        longitude = mass_request.longitude,
                        latitude = mass_request.latitude,
                        minradius = mass_request.minradius,
                        maxradius = mass_request.maxradius)
  receive = get_request(request)
  text = String(receive.body)
  text = strip(text, '#')
  df = CSV.read(IOBuffer(text),delim='|', typemap=Dict(MDateTime=>DateTime,
                MFloat=>Float64,MString=>String))
  newnames =  map(lowercase, String.(names(df)))
  newnames = strip.(newnames)
  newnames = Symbol.(newnames)
  names!(df, newnames)

  return df
end

"""
  mass_download(mass_request, mseed_storage, stationxml_storage)

Launch a data download.

INPUTS TO DO:
  - define domain structure for lat/lon filtering
  - define restrictions structure for time/network/station/channel filtering
"""
function mass_download(mass_request::IRISMassRequest,
                  mseed_storage::String,
                  stationxml_storage::String;
                  download_chunk_size_in_mb::Int=20,
                  threads_per_client::Int=3)
  # function will mirror the obspy mass_downloader.download() function
  # Order of operations TO-DO:
  # 1. Get station availability with station_availability() DONE
  # 2. Setup directory structure based on mseed_storage path DONE
  # 3. Setup directory stucture based on stationxml_storage path DONE
  # 4. Download stationxml for each station - should include response DONE
  # 5. Chunk download into requests based on length of time per mseed file DONE
  #    This should be done in parallel in future using pmap or threads
  # 6. Send HTTP requests to IRIS using IRISTimeSeries DONE
  # 7. Save request to miniseed DONE

  df = station_availability(mass_request)
  # update latest date to today
  latest = DateTime(year(now()),month(now()), day(now()))
  df[:endtime][df[:endtime] .> latest,:] .= latest

  # create directory structure for mseed
  n = size(df,1)
  for ii in 1:n
    current = joinpath(mseed_storage,df[:network][ii],df[:station][ii],
                       df[:channel][ii])
    if isdir(current) == false
      mkpath(current)
    end
  end

  # make directory for stationxml_storage
  if isdir(stationxml_storage) == false
    mkpath(stationxml_storage)
  end

  # download stationxml for each station
  download_stationxml(stationxml_storage,df,starttime=mass_request.starttime,
                      endtime=mass_request.endtime,
                      startbefore=mass_request.startbefore,
                      startafter=mass_request.startafter,
                      endbefore=mass_request.endbefore,
                      endafter=mass_request.endafter,
                      location=mass_request.location)

  # download each network/channel/station
  for ii = 1:n
    NSLC = df[ii,:]
    s, e = NSLC[[:starttime,:endtime]]
    if s < mass_request.starttime; s = mass_request.starttime end
    if e > mass_request.endtime; e = mass_request.endtime end
    station, network, location, channel = NSLC[[:station, :network,
                                                :location, :channel]]
    if ismissing(location); location = "--" end

    # check begin/end dates
    if s != DateTime(year(s),month(s),day(s))
      s = DateTime(year(s),month(s),day(s)) + Dates.Day(1)
    end

    if e != DateTime(year(e),month(e),day(e))
      e = DateTime(year(e),month(e),day(e))
    end

    # start/end date ranges
    start_range = s:Dates.Day(1):e-Dates.Day(1)
    end_range = s+Dates.Day(1) - Dates.Millisecond(1/NSLC[:samplerate][1]
                                                  * 1000) :Dates.Day(1):e

    # download from starttime to endtime in chunks
    for jj=1:length(start_range)
      # get filepath
      filename = get_mseed_name(channel, location, start_range[jj],
                                end_range[jj])
      filepath = joinpath(mseed_storage, network, station, channel, filename)
      download_mseed(station, network, channel,start_range[jj],end_range[jj],
                     filepath)
    end
  end

end

"""
  download_stationxml(args)

Downloads StationXML instrument response.
TODO - verify stationXML are formatted correctly
     - will need to load each file with LightXML
     - this should probably be a separate function
     - make this function threaded with using `threads_per_client`
"""
function download_stationxml(stationxml_storage::String, df::DataFrame;
                             starttime::MDateTime=missing,
                             endtime::MDateTime=missing,
                             startbefore::MDateTime=missing,
                             startafter::MDateTime=missing,
                             endbefore::MDateTime=missing,
                             endafter::MDateTime=missing,
                             location::MString=missing)

  # get all unique network/station combinations
  station_df = unique(df[[:network, :station]])
  n = size(station_df,1)
  for ii in 1:n
    net, sta = station_df[:network][ii], station_df[:station][ii]
    request = FDSNStation(station=sta, network=net,
                        level="response", format="xml", starttime=starttime,
                        endtime=endtime,startbefore=startbefore,
                        startafter=startafter,endbefore=endbefore,
                        endafter=endafter,location=location)
    receive = get_request(request)
    text = String(receive.body)
    stationxml_filename = join([net,sta,"xml"],'.')
    stationxml_path = joinpath(stationxml_storage,stationxml_filename)
    open(stationxml_path,"w") do f
      write(f,text)
    end
    print("Downloaded stationXML for $net.$sta $ii out of $n")
  end
end

"""
   download_mseed(station, network, channel, starttime, endtime, filename)

Download mseed from IRIS.
TO-DO: - try/catch with bad requests
       - print statement with completed requests
       - allow for chunked requests with `download_chunk_size_in_mb`
"""
function download_mseed(station::String, network::String, channel::String,
                      starttime::DateTime,endtime::DateTime,
                      filepath::String)

  try
      request = get_request(IRISTimeSeries(network=network,
                                           station=station,
                                           channel=channel,
                                           location=location,
                                           starttime=starttime,
                                           endtime=endtime,
                                           output="miniseed"))
      text = request.body
      # write to mseed to filepath
      open(filepath, "w") do file
        write(file,text)
      end
      println("Downloaded $network.$station.$location.$channel.$starttime.$endtime")
  catch
      println("Unable to download $network.$station.$location.$channel.$starttime.$endtime")
  end
end

"""
  get_mseed_name(channel, location starttime, endtime)

"""
function get_mseed_name(channel, location, starttime, endtime)
  fmt_string = "%Y%m%dT%H:%M:%S"
  s = Dates.format(starttime, fmt_string)
  e = Dates.format(endtime, fmt_string)
  if location == "--"; location = "" end
  return "$channel.$location.$(starttime)Z.$(endtime)Z.mseed"
  end

network = "CI"
station = "MWC"
channel = "BHZ"
location = "--"
starttime = DateTime(2019, 2, 1, 0, 0, 0)
endtime = now()
req = IRISMassRequest(network=network, station=station, channel=channel,
                      location=location, starttime=starttime,
                      endtime=endtime, output="miniseed")
mseed_storage = "/home/timclements/TEST/MSEED"
stationxml_storage = "/home/timclements/TEST/XML"
mass_download(req,mseed_storage, stationxml_storage)
