# Copied from my old python code, need serious rework.....

function _to_nullable_bool(b::Bool)
    return Nullable(b)
end

function _to_nullable_bool(b::Nullable{Bool})
    return b
end

# There's probably better ways to do this...
function _sprintf(fmt::AbstractString, args...)
    ex = quote
        @sprintf($fmt, $args...)
    end
    return eval(ex)
end

# input: observable, error
# output: formatted observable +- error in scientific notation
function _format_unc(a::Number, s::Number, unit::AbstractString,
                     _sci::Nullable{Bool}, tex::Bool)
    if s <= 0
        return @sprintf("%f%s", a, unit)
    end

    sci = get(_sci, (s >= 100) || max(abs(a), s) < .1)

    la = floor(Int, log10(abs(a)))
    ls = floor(Int, log10(s))
    fs = floor(s * 10.0^(1 - ls))
    if sci
        fa = a * 10.0^(-la)
        dl = la - ls + 1
    else
        fa = a
        dl = 1 - ls
    end
    dl = max(dl, 0)

    ss = if dl == 1
        @sprintf("%.1f", fs / 10)
    else
        @sprintf("%.0f", fs)
    end

    return if sci
        if tex
            _sprintf("%.$(dl)f(%s)\\times10^{%d}{%s}", fa, ss, la, unit)
        else
            _sprintf("%.$(dl)f(%s)*10^%d%s", fa, ss, la, unit)
        end
    else
        _sprintf("%.$(dl)f(%s)%s", fa, ss, unit)
    end
end

function _get_if_list(lst::Union{AbstractString,Number}, idx, _def)
    return lst
end

function _get_if_list(lst::Nullable, idx, _def)
    return _def
end

function _get_if_list(lst, idx, _def)
    try
        return lst[idx]
    catch
        return _def
    end
end

function format_unc(vals::Number, uncs; unit::AbstractString="",
                    sci=Nullable{Bool}(), tex::Bool=false)
    return  _format_unc(vals, uncs::Number, unit, _to_nullable_bool(sci), tex)
end

function format_unc(vals, uncs; unit="", sci=Nullable{Bool}(), tex=false)
    return [_format_unc(v, _get_if_list(uncs, i, 0), _get_if_list(unit, i, ""),
                        _get_if_list(sci, i, Nullable{Bool}()),
                        _get_if_list(tex, i, false))
            for (v, i) in enumerate(vals)]
end
