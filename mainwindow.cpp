#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{   //Initialize Widget
    this->audio_file = nullptr;
    ui->setupUi(this);
    ui->sound_raw->addGraph();
    ui->sound_raw_fft->addGraph();
    ui->sound_filtered->addGraph();
    ui->sound_inversed_filtered->addGraph();
    ui->sound_raw_fft_phase->addGraph();
    ui->sound_phase_filtered->addGraph();
    connect(&this->player, SIGNAL(mediaStatusChanged(QMediaPlayer::MediaStatus)), this, SLOT(on_play()));
}

MainWindow::~MainWindow()
{
    delete ui;
}

    //Function for IFFT
vector<Complex> MainWindow::ifft(const vector<Complex> & data) {
    vector<Complex> result;
    fftw_complex *in, *out;
    in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * data.size());
    for (auto i = 0; i < data.size(); ++i) {
      in[i][0] = data[i].r;
      in[i][1] = data[i].i;
    }
    out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * data.size());
    fftw_plan _plan;
    _plan = fftw_plan_dft_1d(data.size(), in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(_plan);
    for (auto i = 0; i < data.size(); ++i) {
        Complex c;
        c.r = out[i][0];// / (1. * data.size());
        c.i = out[i][1];// / (1. * data.size());
        result.push_back(c);
    }
    fftw_destroy_plan(_plan);
    fftw_free(in);
    fftw_free(out);
    return result;
}
//FFT Function, already used when loading a file.
vector<Complex> MainWindow::fft(const QVector<double> & data) {
  vector<Complex> result;
  fftw_complex *in, *out;
  in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * data.size());
  for (auto i = 0; i < data.size(); ++i) {
    in[i][0] = data[i];
    in[i][1] = .0;
  }
  out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * data.size());
  fftw_plan _plan;
  _plan = fftw_plan_dft_1d(data.size(), in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(_plan);
  for (auto i = 0; i < data.size(); ++i) {
      Complex c;
      c.r = out[i][0];
      c.i = out[i][1];
      result.push_back(c);
  }
  fftw_destroy_plan(_plan);
  fftw_free(in);
  fftw_free(out);
  return result;
}

vector<Complex> MainWindow::lowpass_filter_fft(const double & freq, const vector<Complex> & fft) {
  vector<Complex> result;
  for (auto i = 0; i < fft.size() / 2; i++) {
    double f = i / (fft.size() / (audio_file->fileFormat().sampleRate() * 1.));
    Complex c;
    c.i = 0;
    c.r = 0;
    if (f < freq) {
        c.i = fft[i].i;
        c.r = fft[i].r;
    }
    result.push_back(c);
  }
  for (auto i = fft.size() / 2; i < fft.size(); i++) {
    double f = (fft.size() - i) / (fft.size() / (audio_file->fileFormat().sampleRate() * 1.));
    Complex c;
    c.i = 0;
    c.r = 0;
    if (f < freq) {
        c.i = fft[i].i;
        c.r = fft[i].r;
    }
    result.push_back(c);
  }
  return result;
}

QVector<double> MainWindow::freq_domain(const QVector<double> & data,
                                        double & fmin, double & fmax,
                                        double & vmin, double & vmax) {
    vmin = data[0];
    vmax = data[0];
    QVector<double> items;

    fmin = audio_file->fileFormat().sampleRate();
    fmax = -audio_file->fileFormat().sampleRate();
    qDebug() << "T="<< data.size() / (audio_file->fileFormat().sampleRate() * 1.);

    for (auto i = 0; i < data.size(); i++) {
      if (data[i] < vmin) vmin = data[i];
      if (data[i] > vmax) vmax = data[i];
      double item_val = i / (data.size() / (audio_file->fileFormat().sampleRate() * 1.));

      items.append(item_val);
      if (fmax < item_val) fmax = item_val;
      if (fmin > item_val) fmin = item_val;
    }
    return items;
}
//Play audion function
void MainWindow::play_audio_lr(const QVector<double> & data, const QAudioFormat & format) {
  this->player.stop();
  audio_file->seek(0);
  QByteArray buffer = audio_file->readLine(audio_file->headerLength());
  QDataStream out(&buffer, QIODevice::WriteOnly | QIODevice::Append);
  if (audio_file->fileFormat().byteOrder() == QAudioFormat::LittleEndian)
    out.setByteOrder(QDataStream::LittleEndian);
  else out.setByteOrder(QDataStream::BigEndian);
  out.setByteOrder(QDataStream::BigEndian);
  for (auto d : data) out << static_cast<qint16>(d);
  QBuffer * fw = new QBuffer(this);
  fw->setData(buffer);
  fw->open(QIODevice::ReadOnly);
  this->player.setMedia(QMediaContent(), fw);
}

void MainWindow::plot_amplitude(const QVector<double> & data) {
  double min_f, max_f, min_val, max_val;
  auto items = this->freq_domain(data, min_f, max_f, min_val, max_val);
  ui->sound_raw_fft->xAxis->setRange(min_f, max_f);
  ui->sound_raw_fft->yAxis->setRange(min_val, max_val);
  ui->sound_raw_fft->graph(0)->setData(items, data);
  ui->sound_raw_fft->replot();
}

void MainWindow::plot_amplitude_filtered(const QVector<double> & data) {
    double min_f, max_f, min_val, max_val;
    auto items = this->freq_domain(data, min_f, max_f, min_val, max_val);
    ui->sound_filtered->xAxis->setRange(min_f, max_f);
    ui->sound_filtered->yAxis->setRange(min_val, max_val);
    ui->sound_filtered->graph(0)->setData(items, data);
    ui->sound_filtered->replot();
}

void MainWindow::plot_filtered(const QVector<double> & data) {
    double min_f, max_f, min_val, max_val;
    auto items = this->freq_domain(data, min_f, max_f, min_val, max_val);
    ui->sound_inversed_filtered->xAxis->setRange(min_f, max_f);
    ui->sound_inversed_filtered->yAxis->setRange(min_val, max_val);
    ui->sound_inversed_filtered->graph(0)->setData(items, data);
    ui->sound_inversed_filtered->replot();
}

void MainWindow::plot_phase(const QVector<double> & data) {
  double min_f, max_f, min_val, max_val;
  auto items = this->freq_domain(data, min_f, max_f, min_val, max_val);
  ui->sound_raw_fft_phase->xAxis->setRange(min_f, max_f);
  ui->sound_raw_fft_phase->yAxis->setRange(min_val, max_val);
  ui->sound_raw_fft_phase->graph(0)->setData(items, data);
  ui->sound_raw_fft_phase->replot();
}

void MainWindow::plot_phase_filtered(const QVector<double> & data) {
  double min_f, max_f, min_val, max_val;
  auto items = this->freq_domain(data, min_f, max_f, min_val, max_val);
  ui->sound_phase_filtered->xAxis->setRange(min_f, max_f);
  ui->sound_phase_filtered->yAxis->setRange(min_val, max_val);
  ui->sound_phase_filtered->graph(0)->setData(items, data);
  ui->sound_phase_filtered->replot();
}

void MainWindow::on_load_audio_clicked() {
  audio_filename = QFileDialog::getOpenFileName(this,
    tr("Open Sound"), "/home/", tr("Sound Files (*.wav)"));
  if (audio_filename == "") return;
  if (!audio_file)
    audio_file = new WavFile(this);
  audio_file->open(audio_filename);
  QDataStream data(audio_file);
  if (audio_file->fileFormat().byteOrder() == QAudioFormat::LittleEndian)
    data.setByteOrder(QDataStream::LittleEndian);
  else data.setByteOrder(QDataStream::BigEndian);
  audio_raw.clear();
  audio_amplitude.clear();
  audio_time_stamps.clear();
  audio_fft_angle.clear();
  qint16 min_val=32767, max_val=-32767;
  float rate = audio_file->fileFormat().sampleRate() * 1.;
  while(!data.atEnd()) {
    qint16 x;
    data >> x;
    if (min_val > x) min_val = x;
    if (max_val < x) max_val = x;
    audio_time_stamps.append(audio_raw.size() / rate);
    audio_amplitude.append(static_cast<double>(x));
    audio_raw.append(x);
  }
  ui->sound_raw->xAxis->setRange(0, audio_raw.size() / rate);
  ui->sound_raw->yAxis->setRange(min_val, max_val);
  ui->sound_raw->graph(0)->setData(audio_time_stamps, audio_amplitude);
  ui->sound_raw->replot();
  this->fft_result = this->fft(audio_amplitude);
  for (auto r : this->fft_result) {
    audio_fft_amplitude.append(pow(r.r * r.r + r.i * r.i, .5));

    audio_fft_angle.append(atan2(r.i, r.r) * 180. / M_PI);
  }
  this->plot_amplitude(audio_fft_amplitude);
  this->plot_phase(audio_fft_angle);
}

void MainWindow::on_filter_clicked()
{
  int frequency = ui->freqSelect->value();
  if (this->fft_result.size() == 0) return;
  this->fft_result_filtered = this->lowpass_filter_fft(frequency, this->fft_result);
  audio_fft_amplitude.clear();
  audio_fft_angle.clear();
  for (auto r : this->fft_result_filtered) {
    audio_fft_amplitude.append(pow(r.r * r.r + r.i * r.i, .5));
    audio_fft_angle.append(atan2(r.i, r.r) * 180. / M_PI);
  }
  this->plot_amplitude_filtered(audio_fft_amplitude);
  this->plot_phase_filtered(audio_fft_angle);
}

vector<Complex> MainWindow::highpass_filter_fft(const double & cover, const vector<Complex> & fft) {
    vector<Complex> result;
    QVector2D mid(ui->image->width() / 2., ui->image->height() / 2.);
    for (auto y = 0; y < ui->image->height(); ++y) {
      for (auto x = 0; x < ui->image->width(); ++x) {
          QVector2D current(shift(x, ui->image->width()),
                            shift(y, ui->image->height()));
          Complex c;
          c.i = 0;
          c.r = 0;
          if (mid.distanceToPoint(current) > cover) {
              c.i = fft[x + y * ui->image->height()].i;
              c.r = fft[x + y * ui->image->height()].r;
          }
          result.push_back(c);
      }
    }

    return result;
}

void MainWindow::on_load_image_clicked()
{
   QString image_filename = QFileDialog::getOpenFileName(this,
      tr("Open Image"), "/home/", tr("Image Files (*.bmp *.jpg *.png)"));
   if (image_filename == "") return;
   QPixmap img(image_filename);
   img = img.scaled(ui->image->width(), ui->image->height());
   QImage gs_image(ui->image->width(), ui->image->height(), QImage::Format_RGB32);
   QVector<double> image_raw;

   for (auto y = 0; y < img.height(); ++y) {
     for (auto x = 0; x < img.width(); ++x) {
       QColor col = img.toImage().pixelColor(x, y);
       double val = (col.redF() + col.blueF() + col.greenF()) / 3.;
       QColor gs;
       gs.setRedF(val);
       gs.setBlueF(val);
       gs.setGreenF(val);
       gs_image.setPixelColor(x, y, gs);
       image_raw.push_back(val);
     }
   }
   img = img.fromImage(gs_image);
   ui->image->setPixmap(img);
   this->image_result_fft = this->fft(image_raw);
   this->render_amplitude(this->image_result_fft);
   this->render_phase(this->image_result_fft);
   auto filtered_result = this->highpass_filter_fft(18.0, this->image_result_fft);
   this->render_amplitude_filtered(filtered_result);
   this->render_phase_filtered(filtered_result);
   auto ifft_filtered_result = this->ifft(filtered_result);
   this->render_filtered(ifft_filtered_result);
}

uint32_t MainWindow::shift(const uint32_t & x, const uint32_t & max_val) {
  if (x > max_val / 2) return x - max_val / 2;
  if (x <= max_val / 2) return x + max_val / 2 - 1;
  return x;
}

void MainWindow::on_play() {
    player.play();
}

QImage MainWindow::config_amplitude(const vector<Complex> & data) {
    QImage fft_image(ui->image->width(), ui->image->height(), QImage::Format_RGB32);
    if (data.size() < 1) fft_image;
    auto max_amp = log(pow(data[0].i * data[0].i + data[0].r* data[0].r, 0.5));
    for (auto d : data) {
      auto amp = log(pow(d.i * d.i + d.r * d.r, 0.5));
      if (max_amp < amp) max_amp = amp;
    }
    for (auto y = 0; y < ui->image->height(); ++y) {
      for (auto x = 0; x < ui->image->width(); ++x) {
          QColor gs;
          double r = data[x + y * ui->image->height()].r;
          double i = data[x + y * ui->image->height()].i;
          double amplitude = log(pow(r * r + i * i, 0.5)) / max_amp;
          gs.setRedF(amplitude);
          gs.setBlueF(amplitude);
          gs.setGreenF(amplitude);
          fft_image.setPixelColor(shift(x, ui->image->width()),
                                  shift(y, ui->image->height()), gs);
      }
    }
   return fft_image;
}

void MainWindow::render_amplitude(const vector<Complex> & data) {
    QImage fft_image = this->config_amplitude(data);
    QPixmap img;
    img = img.fromImage(fft_image);
    img = img.scaled(ui->image_fft_amplitude->width(), ui->image_fft_amplitude->height());
    ui->image_fft_amplitude->setPixmap(img);
}

void MainWindow::render_amplitude_filtered(const vector<Complex> &data) {
    QImage fft_image = this->config_amplitude(data);
    QPixmap img;
    img = img.fromImage(fft_image);
    img = img.scaled(ui->image_fft_amplitude_filtered->width(), ui->image_fft_amplitude_filtered->height());
    ui->image_fft_amplitude_filtered->setPixmap(img);
}

void MainWindow::render_filtered(const vector<Complex> & data) {
  QImage fft_image(ui->image->width(), ui->image->height(), QImage::Format_RGB32);
  for (auto y = 0; y < ui->image->height(); ++y) {
    for (auto x = 0; x < ui->image->width(); ++x) {
        QColor gs;
        double r = data[x + y * ui->image->height()].r;
        double amplitude = r / (1. * ui->image->height() * ui->image->width());
        gs.setRedF(amplitude);
        gs.setBlueF(amplitude);
        gs.setGreenF(amplitude);
        fft_image.setPixelColor(x, y, gs);
    }
  }

  QPixmap img;
  img = img.fromImage(fft_image);
  img = img.scaled(ui->image_fft_filtered->width(), ui->image_fft_filtered->height());
  ui->image_fft_filtered->setPixmap(img);
}

QImage MainWindow::config_phase(const vector<Complex> & data) {
    QImage fft_image(ui->image->width(), ui->image->height(), QImage::Format_RGB32);
    auto c = 0;
    for (auto y = 0; y < ui->image->height(); ++y) {
      for (auto x = 0; x < ui->image->width(); ++x) {
          QColor gs;
          double phase_normalized = (atan2(data[c].i, data[c].r) / M_PI + 1.) * 0.5;
          gs.setRedF(phase_normalized);
          gs.setBlueF(phase_normalized);
          gs.setGreenF(phase_normalized);
          fft_image.setPixelColor(shift(x, ui->image->width()),
                                  shift(y, ui->image->height()), gs);
          c++;
      }
    }
    return fft_image;
}

void MainWindow::render_phase(const vector<Complex> & data) {
    QImage fft_image = this->config_phase(data);
    QPixmap img;
    img = img.fromImage(fft_image);
    img = img.scaled(ui->image_fft_phase->width(), ui->image_fft_phase->height());
    ui->image_fft_phase->setPixmap(img);
}

void MainWindow::render_phase_filtered(const vector<Complex> & data) {
    QImage fft_image = this->config_phase(data);
    QPixmap img;
    img = img.fromImage(fft_image);
    img = img.scaled(ui->image_fft_phase_filtered->width(), ui->image_fft_phase_filtered->height());
    ui->image_fft_phase_filtered->setPixmap(img);
}

void MainWindow::on_fft_clicked() {
  if (fft_result_filtered.size() == 0) return;
  auto ifft_results = this->ifft(fft_result_filtered);
  audio_fft_amplitude.clear();
  for (auto r : ifft_results) {
    audio_fft_amplitude.append(r.r / (1. * ifft_results.size()));
  }
  this->plot_filtered(audio_fft_amplitude);
}

void MainWindow::on_play_sound_clicked()
{
    this->play_audio_lr(this->audio_amplitude, this->audio_file->fileFormat());

}

void MainWindow::on_play_sound_filtered_clicked()
{
    this->play_audio_lr(this->audio_fft_amplitude, this->audio_file->fileFormat());
}

void MainWindow::on_Stop_clicked()
{
    this->player.pause();
}
