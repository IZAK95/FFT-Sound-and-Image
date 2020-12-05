#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <fftw3.h>
#include <QFileDialog>
#include <QAudioFormat>
#include <QAudioDecoder>
#include <QDebug>
#include <vector>
#include <QPixmap>
#include <QVector2D>
#include <QAudioBuffer>
#include <QMediaPlayer>
#include <QAudioOutput>
#include "wavfile.h"

using std::vector;
QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE
struct Complex {
  double r;
  double i;
};

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_load_audio_clicked();
    void on_filter_clicked();
    void on_play();
    void on_load_image_clicked();

    void on_fft_clicked();

    void on_play_sound_clicked();

    void on_play_sound_filtered_clicked();

    void on_Stop_clicked();

    void on_spinBox_valueChanged(int arg1);

private:
    vector<Complex> fft(const QVector<double> & data);
    vector<Complex> ifft(const vector<Complex> & data);
    QVector<double> freq_domain(const QVector<double> & data,
                         double & fmin, double &fmax,
                         double &vmin, double & vmax);
    void play_audio_lr(const QVector<double> & data, const QAudioFormat &);
    void plot_amplitude(const QVector<double> & data);
    void plot_phase(const QVector<double> & data);
    void plot_amplitude_filtered(const QVector<double> & data);
    void plot_filtered(const QVector<double> & data);
    void plot_phase_filtered(const QVector<double> & data);
    void render_amplitude(const vector<Complex> & data);
    void render_phase(const vector<Complex> & data);
    void render_amplitude_filtered(const vector<Complex> & data);
    void render_filtered(const vector<Complex> & data);
    void render_phase_filtered(const vector<Complex> & data);
    QImage config_amplitude(const vector<Complex> & data);
    QImage config_phase(const vector<Complex> & data);
    vector<Complex> lowpass_filter_fft(const double & freq, const vector<Complex> & fft);
    vector<Complex> highpass_filter_fft(const double & cover, const vector<Complex> & fft);
    uint32_t shift(const uint32_t & x, const uint32_t & max_val);
    Ui::MainWindow *ui;
    QString audio_filename;
    WavFile * audio_file;
    QVector<qint16> audio_raw;
    vector<Complex> fft_result;
    vector<Complex> fft_result_filtered;
    QVector<double> audio_amplitude;
    QVector<double> audio_time_stamps;
    QVector<double> audio_fft_amplitude;
    QVector<double> audio_fft_angle;
    QVector<double> audio_fft_filtered;
    QByteArray image_data;
    vector<Complex> image_result_fft;
    vector<Complex> image_result_fft_filtered;
    QMediaPlayer player;


};
#endif // MAINWINDOW_H
