import argparse
import numpy as np
import bokeh
from collections import OrderedDict
from bokeh.plotting import figure, output_file, save, show
from bokeh.models import ColumnDataSource, CustomJS, Title, HoverTool, Span, NormalHead, Arrow, LinearColorMapper, ColorBar, NumeralTickFormatter, Range1d, Legend
from bokeh.palettes import viridis, plasma
from bokeh.layouts import column, row, layout
from bokeh.models.widgets import Slider, Div, CheckboxGroup


template = """
{% block preamble %}
<style>
    .acoeffs {
        background-color: #EFF7FE;
        padding: 5px;
        width: inherit;
    }
    .bcoeffs {
        background-color: #E8FFF3;
        padding: 5px;
        width: inherit;
    }
    .tqcoeffs {
        background-color: #EDE1EA;
        padding: 5px;
        width: inherit;
    }
</style>


<script type="text/javascript">
    function linspace(a,b,n) {
        if(typeof n === "undefined") n = Math.max(Math.round(b-a)+1,1);
        if(n<2) { return n===1?[a]:[]; }
        var i,ret = Array(n);
        n--;
        for(i=n;i>=0;i--) { ret[i] = (i*b+(n-i)*a)/n; }
        return ret;
    }

    function polyeval3(coeffs) {
        return function(x){
            return coeffs[2]*x*x + coeffs[1]*x + coeffs[0];
        }
    }

    function polyeval4(coeffs) {
        return function(x){
            return coeffs[3]*x*x*x + coeffs[2]*x*x + coeffs[1]*x + coeffs[0];
        }
    }

    window.normalizeAndCenter=function(ay) {
        let maxval = Math.max(...ay)
        return ay.map(x => x / maxval - 0.5)
    }

    window.splineWave=function(coeffs, time) {
        const wave = time.map(polyeval4(coeffs));
        return wave;
    }
    
    window.counts=function(ay, N) {
        let dx = Math.max(...ay) / N
        ay = ay.map(el => Math.floor(el / dx))
        let counts = new Array(N).fill(0);
        for (var i = 0; i < ay.length; i++) {
            counts[ay[i]] += 1;
        }
        if (counts.length != N) {
            counts.pop();
        }
        return counts;
    }
    
    window.reciprocal=function(ay) {
        return ay.map(x => 1/x);
    }
    
    window.quadWave=function(coeffs, time) {
        const wave = time.map(polyeval3(coeffs));
        return wave;
    }
    window.setSpanPosition=function(span, pos) {
        span.location = pos;
    }
    
    window.integrate=function(time, amp) {
        let moment = [0];
        for (i=0; i < (time.length -1); i++) {
            let dt = time[i+1] - time[i];
            let meanamp = (amp[i+1] - amp[i]) / 2 + amp[i];
            moment.push(dt*meanamp + moment[i]);
        }
        return moment;
    }
</script>
{% endblock %}
"""


def addSpan(fig, name, color='black'):
    span = Span(location=-1,
                dimension='height',
                level='underlay',
                name=name,
                line_color=color,
                line_dash='dashed',
                line_width=2)
    fig.add_layout(span)
    return span


def plotAsym(args):
    p1 = figure(width=600, height=400, toolbar_location=None)
    p2 = figure(width=600, height=200, toolbar_location=None)
    pq = figure(width=600, height=400, toolbar_location=None, x_range=p1.x_range)
    pSamplingDensitySQ = figure(width=600, height=200, toolbar_location=None, title="Sampling density (dwell time)")
    pSamplingDensityST = figure(width=600, height=200, toolbar_location=None)
    pSamplingDensityST.yaxis.axis_label = 'Time per sample [us]'
    pSamplingDensitySQ.yaxis.axis_label = 'Time per sample [us]'
    
    spans = {'echo': addSpan(p1, 'echo', 'chocolate'),
             'mid': addSpan(p1, 'mid', '#6B8E23'),
             'echoM': addSpan(p2, 'echoM', 'chocolate'),
             'midM': addSpan(p2, 'midM', '#6B8E23'),
             'echoQ': addSpan(pq, 'echoQ', 'chocolate'),
             'midQ': addSpan(pq, 'midQ', '#6B8E23')}
    divs = {'acoeffs': Div(css_classes=['acoeffs']), 'bcoeffs': Div(css_classes=['bcoeffs']), 'tqcoeffs': Div(css_classes=['tqcoeffs'])}
    lineCDS = ColumnDataSource({'time': [0, 1], 'amp': [0, 0]})
    trapCDS = ColumnDataSource({'time': [0, 1], 'amp': [0, 0]})
    quadCDS = ColumnDataSource({'time': [0, 1], 'amp': [0, 0]})
    momentCDS = ColumnDataSource({'time': [0, 1], 'quad': [0,1], 'spline': [0,1]})
    densityCDS = ColumnDataSource({'samples':[0, 1], 'trap': [1, 1], 'quad': [1,3], 'spline': [1,2]})
    densityyRange = Range1d(0, 100)
    pSamplingDensityST.y_range = densityyRange
    pSamplingDensitySQ.y_range = densityyRange
    p1.line(x='time', y='amp', source=lineCDS, line_color='orchid', line_width=2)
    p2.line(x='time', y='spline', source=momentCDS, line_color='orchid', line_width=2)
    p2.line(x='time', y='quad', source=momentCDS, line_color='olive', line_width=2)
    pq.line(x='time', y='amp', source=trapCDS, line_color='salmon', line_width=2, line_alpha=.5)
    pq.line(x='time', y='amp', source=quadCDS, line_color='olive', line_width=2)
    r_sdqs = pSamplingDensitySQ.varea_stack(['quad', 'spline'], x='samples', color=("olive", "thistle"), source=densityCDS, legend_label=['Quad', 'Spline'])
    r_sdts = pSamplingDensityST.varea_stack(['trap', 'spline'], x='samples', color=("lightsalmon", "thistle"), source=densityCDS, legend_label=['Trap', 'Spline'])
    pSamplingDensityST.legend.orientation = "horizontal"
    pSamplingDensitySQ.legend.orientation = "horizontal"
    pSamplingDensityST.legend.location = 'top_center'
    pSamplingDensitySQ.legend.location = 'top_center'
    #pSamplingDensityST.line(x='samples', y='quad', color="red", source=densityCDS)
#    pSamplingDensityST.varea_stack(['spline', 'quad'], x='samples', color=("grey", "lightgrey"), source=densityCDS)
    p1.xaxis.axis_label = 'Time [ms]'
    p1.yaxis.axis_label = 'Amplitude'
    p2.yaxis.axis_label = 'k-Space coordinate'
    pq.yaxis.axis_label = 'Amplitude'
    p1.yaxis.formatter = NumeralTickFormatter(format="0.0")
    p2.yaxis.formatter = NumeralTickFormatter(format="0.0")
    pq.yaxis.formatter = NumeralTickFormatter(format="0.0")
    cb = CheckboxGroup(labels=["Mirror positive shifts"], active=[0])
    sliders = OrderedDict({'tas': Slider(start=0, end=2, value=1, step=.1, title="Acquisition start", name='tas'),
                           'ttc': Slider(start=0, end=10, value=3, step=.01, title="Time to center", name='ttc'),
                           'tae': Slider(start=0, end=10, value=5, step=.1, title="Acquisition end", name='tae'),
                           'Mp': Slider(start=0, end=5, value=1, step=1, title="Padding area", name='Mp'),
                           'M': Slider(start=0, end=150, value=75, step=5, title="Sampling area", name='M'),
                           'samples': Slider(start=64, end=256, value=96, step=32, title="Samples", name='samples')})
    sliderCallbackUW = CustomJS(
    args={'figs': {'quad': pq, 'moment': p2, 'spline': p1, 'sdsq': pSamplingDensitySQ, 'sdst': pSamplingDensityST},
          'sliders': sliders,
          'spans': spans,
          'divs': divs,
          'ranges': {'density.y': densityyRange},
          'CDS': {'line': lineCDS, 'trap': trapCDS, 'quad': quadCDS, 'sd': densityCDS, 'moment': momentCDS}, 'checkbox': cb},
    code="""
    let M = sliders['M'].value;
    let Mp = sliders['Mp'].value;
    let tae = sliders['tae'].value;
    let tas = sliders['tas'].value;
    let ttc = sliders['ttc'].value;
    let samples = sliders['samples'].value
    let centerpoint = (tae+tas)/2;
    window.setSpanPosition(spans['echo'], ttc);
    window.setSpanPosition(spans['echoM'], ttc);
    window.setSpanPosition(spans['echoQ'], ttc);
    window.setSpanPosition(spans['mid'], centerpoint);
    window.setSpanPosition(spans['midM'], centerpoint);
    window.setSpanPosition(spans['midQ'], centerpoint);
    let flip = ttc > centerpoint && checkbox.active.length
    if (flip) {
        ttc = centerpoint - (ttc - centerpoint);
    }
    let t1 = ttc - tas;
    let t2 = tae - ttc;
    let g0 = 2 * Mp / tas;
    let Mb1 = g0 * t1;
    let Mb2 = g0 * t2;
    let M1 =  M/2 - Mb1 - Mp;
    let M2 =  M/2 - Mb2 - Mp;
    let a1 = (4*(3*M1*Math.pow(t2,4) - M2*Math.pow(t1,4) + 6*M1*Math.pow(t1,2)*Math.pow(t2,2) - 6*M2*Math.pow(t1,2)*Math.pow(t2,2) + 10*M1*t1*Math.pow(t2,3) - 6*M2*Math.pow(t1,3)*t2))/(Math.pow(t1,2)*Math.pow(t2,2)*(2*t1 + 3*t2)*(t1 + t2))
    let a2 = 12*(M2*t1*t1*t1 + 4*M2*t1*t1*t2 - 4*M1*t1*t2*t2 - M1*t2*t2*t2) / (2*t1*t1*t1*t1*t2*t2 + 3*t1*t1*t1*t2*t2*t2);
    let a3 = 4*(M1*t2*t2*t2*t2 - 2*M2*t1*t1*t1*t1 + t1*t1*t2*t2*(6*M1 - 4*M2) + 5*M1*t1*t2*t2*t2 - 8*M2*t1*t1*t1*t2) / (t1*t1*t1*t1*t2*(t1 + t2)*(3*t2*t2 + 2*t1*t2) );
    let b0 = a3*t1*t1*t1 + a2*t1*t1 + a1*t1;
    let b1 = 3*a3*t1*t1 + 2*a2*t1 + a1;
    let b2 = 3*a3*t1 + a2;
    let b3 = -(3*a3*t1 + a2)/(3*t2);
    let points = 20;
    let times1 = linspace(0, t1, points);
    let times2 = linspace(0, t2, points);
    spline1 = window.splineWave([g0,  a1, a2, a3], times1);
    spline2 = window.splineWave([g0+b0, b1, b2, b3], times2);
    
    let splines = [...spline1, ...spline2]
    if (flip) {
        splines.reverse()
        ttc = sliders['ttc'].value;
    }
    time = [0, tas, ...linspace(tas,ttc,points), ...linspace(ttc, tae, points), tae, tae+tas]
    amp = [0, g0, ...splines, g0, 0]
    CDS['line'].data.amp = amp
    CDS['line'].data.time = time
    CDS['line'].change.emit();

    let trapamp = M/(tae-tas)
    let trapramp = Mp*2/trapamp
    let t_cip = (tae-tas)/2
    let q0 = 1/(2/trapamp - 1/(g0+b0))
    let q2 = 3*M/(2*t_cip*t_cip*t_cip) - 3*(q0)/(t_cip*t_cip)
    CDS['trap'].data.time = [tas-trapramp, tas, tae, tae+trapramp]
    CDS['trap'].data.amp = [0, trapamp, trapamp, 0]
    CDS['trap'].change.emit();

    let timesquad = linspace(-(tae-tas)/2, (tae-tas)/2, points);
    let quad = window.quadWave([q0, 0, q2], timesquad)
    let quadramp = Mp*2/quad[0]
    time = [ tas - quadramp, tas, ...linspace(tas, tae, points), tae + quadramp]
    amp = [ 0, quad[0], ...quad, 0]
    CDS['quad'].data.time = time
    CDS['quad'].data.amp = amp
    CDS['quad'].change.emit()
    
    timesquad = linspace(-(tae-tas)/2, (tae-tas)/2, samples);
    quad = window.quadWave([q0, 0, q2], timesquad)
    quadmoment = window.integrate(timesquad, quad)
    quadmoment = window.normalizeAndCenter(quadmoment)
    
    let trapmoment = window.integrate(linspace(tas, tae, samples), linspace(trapamp, trapamp, samples))
    trapmoment = window.normalizeAndCenter(trapmoment)

    times1 = linspace(0, t1, Math.round( samples*t1/(tae-tas)) );
    times1.pop()
    times2 = linspace(0, t2, Math.round( samples*t2/(tae-tas)) +1 );
    spline1 = window.splineWave([g0,  a1, a2, a3], times1);
    spline2 = window.splineWave([g0+b0, b1, b2, b3], times2);
    let splinemoment = window.integrate(linspace(tas, tae, samples), [...spline1, ...spline2])
    splinemoment = window.normalizeAndCenter(splinemoment)
        
    CDS['moment'].data.spline = splinemoment;
    CDS['moment'].data.quad = quadmoment;
    CDS['moment'].data.time = linspace(tas, tae, samples)
    CDS['moment'].change.emit()
    
    let dwell = 0.0002
    let samplesRaw = Math.round((tae-tas) / dwell)
    console.log(samplesRaw)
    let quadRaw = window.quadWave([q0, 0, q2], linspace(-(tae-tas)/2, (tae-tas)/2, samplesRaw));
    timesquad = linspace(-(tae-tas)/2, (tae-tas)/2, samplesRaw);
    quadmomentRaw = window.integrate(timesquad, quadRaw)
    let quaddwelltime = window.counts(quadmomentRaw, samples)
    
    times1 = linspace(0, t1, Math.round( samplesRaw*t1/(tae-tas)) );
    times1.pop()
    times2 = linspace(0, t2, Math.round( samplesRaw*t2/(tae-tas)) +1 );
    spline1 = window.splineWave([g0,  a1, a2, a3], times1);
    spline2 = window.splineWave([g0+b0, b1, b2, b3], times2);
    let splinemomentRaw = window.integrate(linspace(tas, tae, samplesRaw), [...spline1, ...spline2])
    let splinedwelltime = window.counts(splinemomentRaw, samples)
    
    let trapdwelltime = window.counts(linspace(0, splinemomentRaw.slice(-1)[0], samplesRaw), samples)
    
    quaddwelltime = quaddwelltime.map(x => x * dwell * 1000)
    trapdwelltime = trapdwelltime.map(x => x * dwell * 1000)
    splinedwelltime = splinedwelltime.map(x => x * dwell * 1000)
    CDS['sd'].data.samples = [...Array(quaddwelltime.length).keys()]
    CDS['sd'].data.quad = quaddwelltime
    CDS['sd'].data.trap = trapdwelltime
    CDS['sd'].data.spline = splinedwelltime
    CDS['sd'].change.emit()
    
    ranges['density.y'].start = 0
    ranges['density.y'].end = (Math.max(...trapdwelltime) + Math.max(...splinedwelltime)) * 1.1

    divs['acoeffs'].text = 'g0 = ' + g0.toFixed(2) + '<br>'
                         + 'a1 = ' + a1.toFixed(2) + '<br>'
                         + 'a2 = ' + a2.toFixed(2) + '<br>'
                         + 'a3 = ' + a3.toFixed(2) + '<br>';
    divs['bcoeffs'].text = 's0 =' + (g0+b0).toFixed(2) + '<br>'
                         + 'b1 = ' + b1.toFixed(2) + '<br>'
                         + 'b2 = ' + b2.toFixed(2) + '<br>'
                         + 'b3 = ' + b3.toFixed(2) + '<br>';
    divs['tqcoeffs'].text = 'q0 =' + q0.toFixed(2) + '<br>'
                          + 'q2 =' + q2.toFixed(2) + '<br>'
                          + 'λ = ' + trapamp.toFixed(2) + '<br>'
                          + 'trapdt = ' + (trapdwelltime[samples/2]).toFixed(2) + '<br>'
                          + 'quaddt = ' + (quaddwelltime[samples/2]).toFixed(2) + '<br>'
                          + 'splinedt = ' + (splinedwelltime[samples/2]).toFixed(2) + '<br>';
    sliders['ttc'].title = 'Time to center: (dt = ' + (ttc - centerpoint).toFixed(2) + ' ms)';
    """
    )
    [s.js_on_change('value', sliderCallbackUW) for s in sliders.values()]
    
    
    
    l = layout([
        [p1, [column(list(sliders.values())),
              [cb],
              row(divs['acoeffs'], divs['bcoeffs'], divs['tqcoeffs'])]],
        [pq],
        [p2],
        [pSamplingDensitySQ],
        [pSamplingDensityST]
    ])
    output_file('asym.html')
    save(l, template=template)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--svg',
                        help='SVG output',
                        default=True)
    args = parser.parse_args()
    plotAsym(args)