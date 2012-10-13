#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTimer>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_marker.h>
#include <qwt_symbol.h>
#include <qwt_plot_grid.h>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    
private:
    Ui::MainWindow *ui;
    QwtPlotCurve **curve;
    QwtPlotMarker **mark;
    QTimer *timer;
    QwtPlotGrid *grid;

    double **x, **y, **z;
    double *w;
    double d, a, mu, r, v;
    double Tmin, Tmax, dt;
    int count, num;
    int t;

//    double *vx, *vy;

    void solveStep(int i);

public slots:
    void replot();

signals:
    void stepChanged(int value);

private slots:
    void on_pauseBtn_clicked();
};

#endif // MAINWINDOW_H
