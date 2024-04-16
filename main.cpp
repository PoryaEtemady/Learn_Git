#include "signalfeatures.h"

SignalFeatures::SignalFeatures(QWidget *parent, int channel_number, QColor *plot_colors[])
    :QWidget{parent},m_channel_number(channel_number),m_fft_enable_state(false),m_selected_channel(0),m_fft_powerSpec_state(false)
{
    m_last_p2p_value = 0;
    m_last_average_value =0;
    m_last_max_value = std::numeric_limits<double>::min();;
    m_last_min_value = std::numeric_limits<double>::max();

    for (int cnt_ch = 0; cnt_ch < m_channel_number; ++cnt_ch)
    {
        m_last_fft_max_peak_freq[cnt_ch] = 0;
    }

    connect(this,&SignalFeatures::newPeak2PeakValue,this,&SignalFeatures::stringizeVp2p_FFT);
    connect(this,&SignalFeatures::newFFTValue,this,&SignalFeatures::stringizeVp2p_FFT);

    h_layout = new QHBoxLayout;
    fft_plot   = new QCustomPlot();
    //fft_plot->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

    fft_plot->setInteractions(QCP::iRangeZoom | QCP::iRangeDrag);

    h_layout->addWidget(fft_plot, 1);
    /*
    fft_plot->addGraph();
    QColor color (0, 0, 0);
    fft_plot->graph(0)->setPen(QPen(QColor(color)));
    */
    QColor *color_arr[m_channel_number];
    for(int cnt=0;cnt<m_channel_number;cnt++)
    {
        color_arr[cnt] = new QColor();
        color_arr[cnt] = plot_colors[cnt];
    }

    for(uint8_t cnt=0; cnt<m_channel_number; cnt++)
    {
        fft_plot->addGraph();
        fft_plot->graph(cnt)->setPen(QPen(QColor(*color_arr[cnt])));
    }
    fft_plot->replot(QCustomPlot::rpImmediateRefresh);
    setLayout(h_layout);
    this->setGeometry(0,0,407,320);
    this->setVisible(true);
}

void SignalFeatures::Peak2Peak(const QVector<double> &data)
{
    double max_read = std::numeric_limits<double>::min();
    double min_read = std::numeric_limits<double>::max();
    for (float element : data)
    {
//        if(element<0)
//            element*=-1;
        if(element >max_read)
        {
            max_read = element;
        }
        if(element < min_read)
        {
            min_read = element;
        }
    }

    if(max_read > m_last_max_value)
    {
        m_last_max_value = max_read;
        emit newMaxValue(m_last_max_value);
    }

    if(min_read < m_last_min_value)
    {
        m_last_min_value = min_read;
        emit newMinValue(m_last_min_value);
    }

    m_last_p2p_value = max_read - min_read;
    emit newPeak2PeakValue(m_last_p2p_value);

//    qDebug()<< "V p2p = "<<m_last_p2p_value*1000;
//    qDebug()<< "V max = "<<max_read;
//    qDebug()<< "V min = "<<min_read;
}

void SignalFeatures::Average(const QVector<double> &data)
{

}

void SignalFeatures::Max(const QVector<double> &data)
{
    double max_read = std::numeric_limits<double>::min();
    for (float element : data)
    {
        if(element >max_read)
        {
            max_read = element;
        }
    }
    if(max_read > m_last_max_value)
    {
        m_last_max_value = max_read;
    }
    emit newMaxValue(m_last_max_value);
}

void SignalFeatures::Min(const QVector<double> &data)
{

}

void SignalFeatures::CalcFFT(const QVector<double> &data)
{
    if(m_fft_enable_state)
    {
        //qDebug()<<"Data length = "<<signalLength;
        const int N = 32000;
        double *inputSignal;

// Initialize input signal with some data
#if 0
    static int cnt =0;
    cnt+=100;
    const int signalLength = 32000;
    const int FREQUENCY = cnt;//500;
    inputSignal = (double*) fftw_malloc(sizeof(double) * signalLength);
    for (int i = 0; i < signalLength; i++)
    {
        double t = i/(double)signalLength;
        inputSignal[i] = sin(2 * M_PI * FREQUENCY * t);
        //qDebug()<<i<<" sin("<<t<<")= "<<inputSignal[i];
    }
#else
        const int signalLength = 32000;//data.length();
        inputSignal = (double*) fftw_malloc(sizeof(double) * signalLength);
        for (int i = 0; i < data.length(); i++)
        {
            inputSignal[i] = data[i];//sin(2 * M_PI * FREQUENCY * t);
            //qDebug()<<i<<" sin("<<t<<")= "<<inputSignal[i];
        }
        for (int i = data.length(); i < 32000; i++)
        {

            inputSignal[i] = 0;//sin(2 * M_PI * FREQUENCY * t);
            //qDebug()<<i<<" sin("<<t<<")= "<<inputSignal[i];
        }
#endif
        fftw_complex *outputFFT;
        outputFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1));

        fftw_plan forwardPlan = fftw_plan_dft_r2c_1d(N, inputSignal, outputFFT, FFTW_ESTIMATE);

        fftw_execute(forwardPlan);
        QVector <double>x,y;
        double max = std::numeric_limits<double>::min();
        double amp =0;
        for (int i = 0; i < N/2+1; i++)
        {
            qDebug()<<outputFFT[i][0];
            amp = sqrt(pow(outputFFT[i][0], 2) + pow(outputFFT[i][1], 2));
            //qDebug() << "Frequency: " << i << " Amplitude: " << amp;
            if(amp>max)
            {
                max = amp;
                m_last_fft_max_peak_freq[m_selected_channel] = i;
            }
            x<<i;
            y<<amp;
        }

        //qDebug()<<"fft_max = "<<max<<" at : "<<m_last_fft_max_peak_freq;

        emit newFFTValue(m_last_fft_max_peak_freq[m_selected_channel]);

        fftw_destroy_plan(forwardPlan);
        fftw_free(inputSignal);
        fftw_free(outputFFT);

#if 1
        fft_plot->xAxis->setRange(0,16000);
        fft_plot->yAxis->setRange(-0.1,max*1.01);//2500
        fft_plot->xAxis->setLabel("Frequency[Hz]");
        fft_plot->yAxis->setLabel("|Amp|");
        fft_plot->graph(0)->setData(x, y);
        fft_plot->replot(QCustomPlot::rpImmediateRefresh);
#endif
    }
}

void SignalFeatures::CalcFFTS(const QVector<double> data[])
{

    if(m_fft_enable_state)
    {
        //qDebug()<<data->
        const int N = 32000;
        double *inputSignal[m_channel_number];
        // Initialize input signal with some data
        const int signalLength = 32000;//data.length();
        fftw_complex *outputFFT[m_channel_number];
        static bool write_once = false;

        if(!write_once)
        {
            SaveDataCSV save_data_csv("1signals",data[0]/*,true*/);
            write_once = true;
        }

        for(int ch_cnt=0;ch_cnt<m_channel_number;ch_cnt++)
        {
            if(m_channels_visibility_state[ch_cnt])
            {
                inputSignal[ch_cnt] = (double*) fftw_malloc(sizeof(double) * signalLength);
                outputFFT[ch_cnt] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1));
                for (int data_index = 0; data_index < data[ch_cnt].length(); data_index++)
                {
                    inputSignal[ch_cnt][data_index] = data[ch_cnt].at(data_index);//sin(2 * M_PI * FREQUENCY * t);
                    //qDebug()<<i<<" sin("<<t<<")= "<<inputSignal[i];
                }
                for (int data_index = data[ch_cnt].length(); data_index < 32000; data_index++)
                {

                    inputSignal[ch_cnt][data_index] = 0;//sin(2 * M_PI * FREQUENCY * t);
                    //qDebug()<<i<<" sin("<<t<<")= "<<inputSignal[i];
                }
                fftw_plan forwardPlan = fftw_plan_dft_r2c_1d(N, inputSignal[ch_cnt], outputFFT[ch_cnt], FFTW_ESTIMATE);
                fftw_execute(forwardPlan);
            }
        }



        QVector <double>x[m_channel_number],y[m_channel_number];
        double max[m_channel_number];
        double amp[m_channel_number];
        static double maxOfAll =std::numeric_limits<double>::min();

        for(int cnt=0;cnt<m_channel_number;cnt++)
        {
            max[cnt] = std::numeric_limits<double>::min();
            amp[cnt] = 0;
        }


        static bool first_cpy_arr[8] = {true,true,true,true,true,true,true,true};
        for(int ch_cnt=0;ch_cnt<m_channel_number;ch_cnt++)
        {
            if(m_channels_visibility_state[ch_cnt])
            {
                for (int i = 0; i < N/2+1; i++)
                {
                    //qDebug()<<"outputFFT ["<<i<<"]= "<<outputFFT[2][i][0];
                    amp[ch_cnt] = sqrt(pow(outputFFT[ch_cnt][i][0], 2) + pow(outputFFT[ch_cnt][i][1], 2));
                    //qDebug() << "Frequency: " << i << " Amplitude: " << amp;
                    if(amp[ch_cnt]>max[ch_cnt])
                    {
                        max[ch_cnt] = amp[ch_cnt];
                        m_last_fft_max_peak_freq[ch_cnt] = i;
                        if(ch_cnt == 0)
                            qDebug()<<"max = "<<max[ch_cnt]/*/data[ch_cnt].length()*/<<"m_last_fft_max_peak_freq[ch_cnt] = "<<i<<" .";
                    }
                    //if(ch_cnt == m_selected_channel)
                    //{
                        if(first_cpy_arr[ch_cnt])
                        {
                            amp_power_spectrum[ch_cnt]<<amp[ch_cnt];
                        }
                        else
                        {
                            if(amp[ch_cnt]>amp_power_spectrum[ch_cnt].at(i))
                            {
                                //qDebug()<<"Before freq ["<<i<<"] = "<<amp_power_spectrum[i];
                                amp_power_spectrum[ch_cnt][i] = amp[ch_cnt];
                                //qDebug()<<amp[ch_cnt];
                                //qDebug()<<"After freq  ["<<i<<"] = "<<amp_power_spectrum[i];
                            }
                        }
                    //}
                    x[ch_cnt]<<i;
                    y[ch_cnt]<<amp[ch_cnt];
                }
                first_cpy_arr[ch_cnt] = false;
            }
        }

        //qDebug()<<"fft_max = "<<max<<" at : "<<m_last_fft_max_peak_freq;
        //emit newFFTValue(m_last_fft_max_peak_freq);
        //fftw_destroy_plan(forwardPlan);
        for (int ch_cnt = 0; ch_cnt < m_channel_number; ch_cnt++)
        {
            if(max[ch_cnt]>maxOfAll)
                maxOfAll = max[ch_cnt];
        }
        fft_plot->xAxis->setRange(0,16000);
        fft_plot->yAxis->setRange(-0.1,maxOfAll*1.01);//2500
        //qDebug()<<"MAxOFAll = "<<maxOfAll;
        fft_plot->xAxis->setLabel("Frequency[Hz]");
        fft_plot->yAxis->setLabel("|Amp|");
        for(int ch_cnt=0;ch_cnt<m_channel_number;ch_cnt++)
        {
            if(m_channels_visibility_state[ch_cnt])
            {
                static int last_selected_ch = m_selected_channel;
                static QCPLayer *default_layer = fft_plot->graph(0)->layer();
                if(!m_fft_powerSpec_state)//ch_cnt!=m_selected_channel) // for FFT Plot
                {
                    fft_plot->graph(ch_cnt)->setVisible(true);
                    fft_plot->graph(ch_cnt)->setData(x[ch_cnt], y[ch_cnt]);
                    fft_plot->replot(QCustomPlot::rpImmediateRefresh);
                    fftw_free(inputSignal[ch_cnt]);
                    fftw_free(outputFFT[ch_cnt]);
                }
                else // for PowerSpectrum Plot
                {
                    static double selected_ch_max[8];
                    if(max[ch_cnt]>selected_ch_max[ch_cnt])
                        selected_ch_max[ch_cnt] = max[ch_cnt];
                    //selected_ch_max[ch_cnt]
                    //fft_plot->graph(ch_cnt)->setLayer()
                    fft_plot->graph(ch_cnt)->setVisible(true);
                    fft_plot->graph(ch_cnt)->setData(x[ch_cnt], amp_power_spectrum[ch_cnt]);
                    fft_plot->graph(m_selected_channel)->setLayer("overlay");
                    if(last_selected_ch != m_selected_channel)
                    {
                        fft_plot->graph(last_selected_ch)->setLayer(default_layer);
                        last_selected_ch = m_selected_channel;
                    }


                    fft_plot->replot(QCustomPlot::rpImmediateRefresh);
                    fftw_free(inputSignal[ch_cnt]);
                    fftw_free(outputFFT[ch_cnt]);
                }
            }
            else
            {
                fft_plot->graph(ch_cnt)->setVisible(false);
            }
        }
    }

}

void SignalFeatures::channelsVisibilityState(const bool visibilityState[])
{
    for (int cnt = 0; cnt < m_channel_number; cnt++)
    {
        m_channels_visibility_state[cnt] = visibilityState[cnt];
        //qDebug()<<"channelsVisibilityState ["<<cnt<<"] Slot : "<<m_channels_visibility_state[cnt];
    }
}

void SignalFeatures::fftEnableState(bool enable)
{
    m_fft_enable_state = enable;
}

float SignalFeatures::last_min_value() const
{
    return m_last_min_value;
}

void SignalFeatures::setLast_min_value(float newLast_min_value)
{
    m_last_min_value = newLast_min_value;
}

void SignalFeatures::stringizeVp2p_FFT(void)
{
    QString strVp2p_FFT="";
    QString str_vp2p="";
    QString str_fft="";

    if(m_last_p2p_value>1)
        str_vp2p = "  Vp_p = " + QString::number(m_last_p2p_value)+" V";
    else if(m_last_p2p_value<1 && m_last_p2p_value >0.001)
        str_vp2p = "  Vp_p = " + QString::number(m_last_p2p_value*1000)+" mV";
    else
        str_vp2p = "  Vp_p = " + QString::number(m_last_p2p_value*1000000)+" uV";

    if(m_fft_enable_state)
    {
        float last_fft_max_peak_freq_selected_channel = m_last_fft_max_peak_freq[m_selected_channel];
        if(last_fft_max_peak_freq_selected_channel<1000)
            str_fft = "  Freq = " + QString::number(last_fft_max_peak_freq_selected_channel) + " Hz";
        else if((last_fft_max_peak_freq_selected_channel>=1000) && (last_fft_max_peak_freq_selected_channel<1000000))
            str_fft = "  Freq = " + QString::number(last_fft_max_peak_freq_selected_channel/1000) + " KHz";
        else if((last_fft_max_peak_freq_selected_channel>=1000000) && (last_fft_max_peak_freq_selected_channel<1000000000))
            str_fft = "  Freq = " + QString::number(last_fft_max_peak_freq_selected_channel/1000) + " MHz";
    }
    strVp2p_FFT = str_vp2p + str_fft;
    emit newVP2P_FFT(strVp2p_FFT);
}

void SignalFeatures::selectedChannelSetter(int selected_channel)
{
    m_selected_channel = selected_channel;
}

void SignalFeatures::setFFT_PowerSpec(bool state)
{
    m_fft_powerSpec_state = state;
}

float SignalFeatures::last_max_value() const
{
    return m_last_max_value;
}

void SignalFeatures::setLast_max_value(float newLast_max_value)
{
    m_last_max_value = newLast_max_value;
}

float SignalFeatures::last_average_value() const
{
    return m_last_average_value;
}

void SignalFeatures::setLast_average_value(float newLast_average_value)
{
    m_last_average_value = newLast_average_value;
}

float SignalFeatures::last_p2p_value() const
{
    return m_last_p2p_value;
}

void SignalFeatures::setLast_p2p_value(float newLast_p2p_value)
{
    m_last_p2p_value = newLast_p2p_value;
}
