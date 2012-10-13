#include "mainwindow.h"
#include "ui_mainwindow.h"

#define F_WIDTH     100
#define F_HEIGHT    100
#define PI          3.14159265

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    num = 10;
    curve = new QwtPlotCurve * [num];
    mark = new QwtPlotMarker * [num];

    for (int j = 0; j < num; ++j)
    {
        curve[j] = new QwtPlotCurve;
        curve[j]->setRenderHint(QwtPlotItem::RenderAntialiased);

        mark[j] = new QwtPlotMarker;
        mark[j]->setSymbol(new QwtSymbol(QwtSymbol::Ellipse, QColor(Qt::transparent), QColor(Qt::black), QSize(5,5)));
    }

    ui->qwtPlot->setAxisScale(QwtPlot::xBottom, 0, F_WIDTH);
    ui->qwtPlot->setAxisScale(QwtPlot::yLeft, 0, F_HEIGHT);
    grid = new QwtPlotGrid;
    grid->setMajPen(QPen(Qt::black,0,Qt::DotLine));
    grid->attach(ui->qwtPlot);

    v = 2;
    d = 0;//0.5;
    Tmin = 0;
    Tmax = 200;
    count = 100000;
    t = 0;
    connect(this, SIGNAL(stepChanged(int)), ui->lcdNumber, SLOT(display(int)));
    ui->lcdNumber->setNumDigits(trunc(log10(count)));

    x = new double * [num];
    y = new double * [num];

//    vx = new double [num];
//    vy = new double [num];

    for (int j = 0; j < num; ++j)
    {
        x[j] = new double[count];
        y[j] = new double[count];

        x[j][0] = (double)rand()/RAND_MAX * F_WIDTH;
        y[j][0] = (double)rand()/RAND_MAX * F_HEIGHT;
    }

    timer = new QTimer(this);
    connect(timer, SIGNAL(timeout()), this, SLOT(replot()));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::solveStep(int i)
{
    for (int j = 0; j < num; ++j)
    {
        double theta = (double)rand()/RAND_MAX * 2 * PI - PI;
        double vx = v * cos(theta);
        double vy = v * sin(theta);

        x[j][i] = x[j][i-1] + vx;
        y[j][i] = y[j][i-1] + vy;

        if (x[j][i] < 0)
        {
            x[j][i] = -x[j][i];
            vx = -vx;
        }

        if (x[j][i] > F_WIDTH)
        {
            x[j][i] = 2 * F_WIDTH - x[j][i];
            vx = -vx;
        }

        if (y[j][i] < 0)
        {
            y[j][i] = -y[j][i];
            vy = -vy;
        }

        if (y[j][i] > F_HEIGHT)
        {
            y[j][i] = 2 * F_HEIGHT - y[j][i];
            vy = -vy;
        }
    }
}

void MainWindow::replot()
{
    emit stepChanged(t);
    t++;
    solveStep(t);

    for (int j = 0; j < num; ++j)
    {
        QVector <double> xr, yr;
        int ind = std::max(0, t - 3000 / (int)Tmax);
        for (int l = ind; l < t; l++)
        {
            xr.append(x[j][l]);
            yr.append(y[j][l]);
        }

        curve[j]->setSamples(xr, yr);
        curve[j]->attach(ui->qwtPlot);
        mark[j]->setValue(x[j][t-1], y[j][t-1]);
        mark[j]->attach(ui->qwtPlot);
    }
    ui->qwtPlot->replot();

    if (t > count - 2)
    {
        killTimer(timer->timerId());
        ui->pauseBtn->setEnabled(false);
        return;
    }
}

void MainWindow::on_pauseBtn_clicked()
{
    if (ui->pauseBtn->text() == "Pause")
    {
        timer->stop();
        ui->pauseBtn->setText("Play");
    }
    else
    {
        timer->start(50);
        ui->pauseBtn->setText("Pause");
    }
}
