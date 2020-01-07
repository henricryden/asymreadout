import argparse
import numpy as np
import bokeh
from collections import OrderedDict
from bokeh.plotting import figure, output_file, save, show
from bokeh.models import ColumnDataSource, CustomJS, Title, HoverTool, Span, NormalHead, Arrow, LinearColorMapper, ColorBar, NumeralTickFormatter
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

    function polyeval4(coeffs) {
        return function(x){
            return coeffs[3]*x*x*x + coeffs[2]*x*x + coeffs[1]*x + coeffs[0];
        }
    }

    window.splineWave=function(coeffs, time) {
        const map1 = time.map(polyeval4(coeffs));
        return map1;
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
    p1 = figure(width=600, height=600, toolbar_location=None)
    p2 = figure(width=600, height=200, toolbar_location=None)
    spans = {'echo': addSpan(p1, 'echo', 'chocolate'),
             'mid': addSpan(p1, 'mid', '#6B8E23'),
             'echoM': addSpan(p2, 'echoM', 'chocolate'),
             'midM': addSpan(p2, 'midM', '#6B8E23')}
    divs = {'acoeffs': Div(css_classes=['acoeffs']), 'bcoeffs': Div(css_classes=['bcoeffs'])}
    lineCDS = ColumnDataSource({'time': [0, 1], 'amp': [0, 0], 'moment': [0, 0]})
    p1.line(x='time', y='amp', source=lineCDS, line_color='navy', line_width=1.5)
    p2.line(x='time', y='moment', source=lineCDS, line_color='navy', line_width=1.5)
    p1.xaxis.axis_label = 'Time [ms]'
    p1.yaxis.axis_label = 'Amplitude'
    p2.yaxis.axis_label = 'k-Space coordinate'
    p1.yaxis.formatter = NumeralTickFormatter(format="0.0")
    p2.yaxis.formatter = NumeralTickFormatter(format="0.0")
    cb = CheckboxGroup(labels=["Mirror positive shifts"], active=[0])
    sliders = OrderedDict({'tas': Slider(start=0, end=2, value=1, step=.1, title="Acquisition start", name='tas'),
                           'ttc': Slider(start=0, end=10, value=3, step=.01, title="Time to center", name='ttc'),
                           'tae': Slider(start=0, end=10, value=5, step=.1, title="Acquisition end", name='tae'),
                           'Mp': Slider(start=0, end=5, value=2, step=.1, title="Padding area", name='Mp'),
                           'M': Slider(start=0, end=50, value=25, step=.1, title="Sampling area", name='M')})
    sliderCallbackUW = CustomJS(
    args={'sliders': sliders, 'spans': spans, 'divs': divs, 'CDS': {'line': lineCDS}, 'checkbox': cb},
    code="""
    let M = sliders['M'].value;
    let Mp = sliders['Mp'].value;
    let tae = sliders['tae'].value;
    let tas = sliders['tas'].value;
    let ttc = sliders['ttc'].value;
    let centerpoint = (tae+tas)/2;
    window.setSpanPosition(spans['echo'], ttc);
    window.setSpanPosition(spans['echoM'], ttc);
    window.setSpanPosition(spans['mid'], centerpoint);
    window.setSpanPosition(spans['midM'], centerpoint);
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
    let moment = window.integrate(time, amp);
    let momentMax = Math.max(...moment);
    moment = moment.map(x => x / momentMax - 0.5);
    CDS['line'].data.moment = moment;
    CDS['line'].change.emit();
    
    divs['acoeffs'].text = 'g0 = ' + g0.toFixed(2) + '<br>'
                         + 'a1 = ' + a1.toFixed(2) + '<br>'
                         + 'a2 = ' + a2.toFixed(2) + '<br>'
                         + 'a3 = ' + a3.toFixed(2) + '<br>';
    divs['bcoeffs'].text = 'b0 =' + b0.toFixed(2) + '<br>'
                         + 'b1 = ' + b1.toFixed(2) + '<br>'
                         + 'b2 = ' + b2.toFixed(2) + '<br>'
                         + 'b3 = ' + b3.toFixed(2) + '<br>';
    sliders['ttc'].title = 'Time to center: (dt = ' + (ttc - centerpoint).toFixed(2) + ' ms)';
    """
    )
    [s.js_on_change('value', sliderCallbackUW) for s in sliders.values()]
    
    
    
    l = layout([
        [p1, [column(list(sliders.values())),
              [cb],
              row(divs['acoeffs'], divs['bcoeffs'])]],
        [p2]
    ])
    save(l, filename='asym.html', template=template)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--svg',
                        help='SVG output',
                        default=True)
    args = parser.parse_args()
    plotAsym(args)