#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cstdlib>
#include <stdio.h>
#include <qcustomplot.h>

using namespace std;

double M1, M2, R1, R2, L1, L2, R, q, X, Y, X1, X2, Y1, Y2, Z, Z1, Z2, inc, A1, A2, dA1, K, dA2, Phi, Theta1, Theta2, p, residu, residu2;
QString QPath;
string path, line, eins, zwei, drei;
QVector <double> phim(1), mag(1), element(8);

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    this->setWindowTitle("LightCurve");

    ui->doubleSpinBox->setValue(0.6);
    ui->doubleSpinBox_2->setValue(0.1);
    ui->doubleSpinBox_3->setValue(0.4);
    ui->doubleSpinBox_4->setValue(0.3);
    ui->doubleSpinBox_5->setValue(0.6);
    ui->doubleSpinBox_6->setValue(0.2);
    ui->doubleSpinBox_7->setValue(1.0);
    ui->doubleSpinBox_8->setValue(85);
    ui->doubleSpinBox_9->setValue(5.0);

    QFont legendFont = font();
    legendFont.setPointSize(25);
    ui->plot->xAxis->setLabelFont(legendFont);
    ui->plot->yAxis->setLabelFont(legendFont);
    ui->plot->xAxis->setTickLabelFont(legendFont);
    ui->plot->yAxis->setTickLabelFont(legendFont);
    ui->plot->axisRect()->setupFullAxesBox(true);
    ui->plot->yAxis->setRange(0,1.2);
    ui->plot->xAxis->setRange(0,1.5);
    ui->plot->xAxis->setLabel("phase");
    ui->plot->yAxis->setLabel("normalized flux");

    ui->lineEdit->setText("phases.dat");
    ui->lineEdit_2->setText("light.dat");
    ui->lineEdit_3->setText("/home/daniels/LightCurve");
    ui->lineEdit_4->setText("messung.dat");


}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_pushButton_2_clicked()
{
    ui->plot->clearGraphs();
    QPath = ui->lineEdit_3->text();
    path = QPath.toUtf8().constData();

    M1 = ui->doubleSpinBox->value();
    M2 = ui->doubleSpinBox_2->value();
    R1 = ui->doubleSpinBox_3->value();
    R2 = ui->doubleSpinBox_4->value();
    L1 = ui->doubleSpinBox_5->value();
    L2 = ui->doubleSpinBox_6->value();
    R = ui->doubleSpinBox_7->value();
    inc = ui->doubleSpinBox_8->value()*M_PI/180;
    K = ui->doubleSpinBox_9->value();

    q=M2/M1;

    ui->doubleSpinBox_10->setValue(q);

    QString input=ui->lineEdit->text();
    string data = input.toUtf8().constData();
    std::ostringstream datNameStream(data);
    datNameStream<<path<<"/"<<data;
    std::string datName = datNameStream.str();
    ifstream dat(datName.c_str());

    int lines=0;

    while(std::getline(dat, line))
       ++ lines;

    dat.clear();
    dat.seekg(0, ios::beg);

    QVector<double> phi(lines);

    for (int i =0; i < lines; i++){
    dat >> eins;
    istringstream istr(eins);
    istr >> phi[i];

    }
    dat.close();

    QVector<double> I(lines);

    QString output=ui->lineEdit_2->text();
    string out = output.toUtf8().constData();
    std::ostringstream outNameStream(out);
    outNameStream<<path<<"/"<<out;
    std::string outName = outNameStream.str();
    ofstream odat(outName.c_str());

    for(int i =0; i<lines; i++){

        X=R*sin(2*M_PI*phi[i]);
        X1=-X/(1.0+(1.0/q));
        X2=X/(1.0+q);
        Y=R*cos(inc)*cos(2*M_PI*phi[i]);
        Y1=-Y/(1+(1/q));
        Y2=Y/(1+q);
        Z=R*sin(inc)*cos(2*M_PI*phi[i]);
        Z1=-Z/(1.0+(1.0/q));
        Z2=Z/(1.0+q);
        p=sqrt(pow((X2-X1),2)+pow((Y2-Y1),2));
        Theta1=2*acos((R2*R2-R1*R1-p*p)/(-2*R1*p));
        Theta2=2*acos((R1*R1-R2*R2-p*p)/(-2*R2*p));

        if(Theta1>M_PI){
            Theta1 = 2*M_PI-Theta1;
        }

        if(Theta2>M_PI){
            Theta2 = 2*M_PI-Theta2;
        }

        dA1=R1*R1*0.5*(Theta1-sin(Theta1));
        dA2=R2*R2*0.5*(Theta2-sin(Theta2));

        if(Z1>Z2){
            if(p>R1+R2){
                A1 = R1*R1*M_PI;
                A2 = R2*R2*M_PI;
            }
            if((p<R1+R2) & (p>sqrt(R1*R1-R2*R2))){
                A1 = R1*R1*M_PI;
                A2 = R2*R2*M_PI-dA1-dA2;
            }
            if((p<sqrt(R1*R1-R2*R2)) & (p>R1-R2)){
                A1 = R1*R1*M_PI;
                A2 = dA2-dA1;
            }
            if(p<=R1-R2){
                A1 = R1*R1*M_PI;
                A2 = 0;
            }
        }
        if(Z1<Z2){
            if(p>R1+R2){
                A1 = R1*R1*M_PI;
                A2 = R2*R2*M_PI;
            }
            if((p<R1+R2) & (p>sqrt(R1*R1-R2*R2))){
                A1 = R1*R1*M_PI-dA1-dA2;
                A2 = R2*R2*M_PI;
            }
            if((p<sqrt(R1*R1-R2*R2)) & (p>R1-R2)){
                A1 = R1*R1*M_PI-R2*R2*M_PI+dA2-dA1;
                A2 = R2*R2*M_PI;
            }
            if(p<=R1-R2){
               A1 = R1*R1*M_PI-R2*R2*M_PI;
               A2 = R2*R2*M_PI;
            }
        }

        I[i] = K*(L1*A1/(R1*R1)+L2*A2/(R2*R2))/(4*M_PI);

        odat<<phi[i]<<"\t"<<I[i]<<endl;
    }
    odat.close();

    I.resize(100);
    phi.resize(100);

    QVector<double> I1(100), I2(100), Ratio(100);

    for(int i =0; i<100; i++){
        phi[i]=i*1.5/100;

        X=R*sin(2*M_PI*phi[i]);
        X1=-X/(1.0+(1.0/q));
        X2=X/(1.0+q);
        Y=R*cos(inc)*cos(2*M_PI*phi[i]);
        Y1=-Y/(1+(1/q));
        Y2=Y/(1+q);
        p=sqrt(pow((X2-X1),2)+pow((Y2-Y1),2));
        Z=R*sin(inc)*cos(2*M_PI*phi[i]);
        Z1=-Z/(1.0+(1.0/q));
        Z2=Z/(1.0+q);
        Theta1=2*acos((R2*R2-R1*R1-p*p)/(-2*R1*p));
        Theta2=2*acos((R1*R1-R2*R2-p*p)/(-2*R2*p));

        if(Theta1>M_PI){
            Theta1 = 2*M_PI-Theta1;
        }

        if(Theta2>M_PI){
            Theta2 = 2*M_PI-Theta2;
        }

        dA1=R1*R1*0.5*(Theta1-sin(Theta1));
        dA2=R2*R2*0.5*(Theta2-sin(Theta2));

        cout<<Theta1<<"\t"<<Theta2<<"\t";

        if(Z1>Z2){
            if(p>R1+R2){
                A1 = R1*R1*M_PI;
                A2 = R2*R2*M_PI;
                cout<<"A1"<<"\t";
            }
            if((p<R1+R2) & (p>sqrt(R1*R1-R2*R2))){
                A1 = R1*R1*M_PI;
                A2 = R2*R2*M_PI-dA1-dA2;
                cout<<"A2"<<"\t";
            }
            if((p<sqrt(R1*R1-R2*R2)) & (p>R1-R2)){
                A1 = R1*R1*M_PI;
                A2 = dA2-dA1;
                cout<<"A3"<<"\t";
            }
            if(p<=R1-R2){
                A1 = R1*R1*M_PI;
                A2 = 0;
                cout<<"A4"<<"\t";
            }
        }
        if(Z1<Z2){
            if(p>R1+R2){
                A1 = R1*R1*M_PI;
                A2 = R2*R2*M_PI;
                cout<<"B1"<<"\t";
            }
            if((p<R1+R2) & (p>sqrt(R1*R1-R2*R2))){
                A1 = R1*R1*M_PI-dA1-dA2;
                A2 = R2*R2*M_PI;
                cout<<"B2"<<"\t";
            }
            if((p<sqrt(R1*R1-R2*R2)) & (p>R1-R2)){
                A1 = R1*R1*M_PI-R2*R2*M_PI+dA2-dA1;
                A2 = R2*R2*M_PI;
                cout<<"B3"<<"\t";
            }
            if(p<=R1-R2){
               A1 = R1*R1*M_PI-R2*R2*M_PI;
               A2 = R2*R2*M_PI;
               cout<<"B4"<<"\t";
            }
        }

        I1[i] = K*L1*A1/(R1*R1)/(4*M_PI);
        I2[i] = K*L2*A2/(R2*R2)/(4*M_PI);
        Ratio[i] = I2[i]/I1[i];
        I[i] = K*(L1*A1/(R1*R1)+L2*A2/(R2*R2))/(4*M_PI);
        cout<<A1<<"\t"<<A2<<"\t"<<I[i]<<endl;

    }

    QString input2=ui->lineEdit_4->text();
    string data2 = input2.toUtf8().constData();
    std::ostringstream dat2NameStream(data2);
    dat2NameStream<<path<<"/"<<data2;
    std::string dat2Name = dat2NameStream.str();
    ifstream dat2(dat2Name.c_str());

    lines=0;

    while(std::getline(dat2, line))
       ++ lines;

    dat2.clear();
    dat2.seekg(0, ios::beg);

    phim.resize(lines);
    mag.resize(lines);

    for (int i =0; i < lines; i++){
    dat2 >> eins >> zwei;
    istringstream istr3(eins);
    istr3 >> phim[i];
    istringstream istr4(zwei);
    istr4 >> mag[i];
    //istringstream istr3(drei);
    //istr3 >> error[i];
    }
    dat2.close();

    QPen pen;
    pen.setWidth(2);
    pen.setColor(Qt::red);
    QPen pen2;
    pen2.setWidth(2);
    pen2.setColor(Qt::blue);
    QPen pen3;
    pen3.setWidth(2);
    pen3.setColor(Qt::green);
    QPen pen4;
    pen4.setWidth(2);
    pen4.setColor(Qt::black);
    QPen pen5;
    pen5.setWidth(2);
    pen5.setColor(Qt::magenta);
    ui->plot->legend->setVisible(true);
    ui->plot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignBottom|Qt::AlignRight);
    ui->plot->addGraph();
    ui->plot->graph(0)->setData(phi, I);
    ui->plot->yAxis->rescale();
    ui->plot->graph(0)->setName("I");
    ui->plot->addGraph();
    ui->plot->graph(1)->setData(phi, I1);
    ui->plot->yAxis->rescale(true);
    ui->plot->graph(1)->setName("I1");
    ui->plot->addGraph();
    ui->plot->graph(2)->setData(phi, I2);
    ui->plot->yAxis->rescale(true);
    ui->plot->graph(2)->setName("I2");
    ui->plot->addGraph();
    ui->plot->graph(3)->setData(phi, Ratio);
    ui->plot->yAxis->rescale(true);
    ui->plot->graph(3)->setName("Ratio");
    ui->plot->addGraph();
    ui->plot->graph(4)->setLineStyle(QCPGraph::lsNone);
    ui->plot->graph(4)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 2));
    ui->plot->graph(4)->addData(phim, mag);
    ui->plot->yAxis->rescale(true);
    ui->plot->graph(4)->setName("Measurement");

   // QCPErrorBars *errorBars = new QCPErrorBars(ui->plot->yAxis);
   // errorBars->removeFromLegend();
   // errorBars->setDataPlottable(ui->plot->graph(4));

    ui->plot->graph(0)->setPen(pen);
    ui->plot->graph(1)->setPen(pen2);
    ui->plot->graph(2)->setPen(pen3);
    ui->plot->graph(3)->setPen(pen4);
    ui->plot->graph(4)->setPen(pen5);

    ui->plot->replot();


}

//*************************************
// fitting lightcurve
//*************************************
void MainWindow::on_pushButton_3_clicked()
{
    QPath = ui->lineEdit_3->text();
    path = QPath.toUtf8().constData();

    residu2=-1;

    M1 = ui->doubleSpinBox->value();
    M2 = ui->doubleSpinBox_2->value();
    element[1] = ui->doubleSpinBox_3->value();
    element[2] = ui->doubleSpinBox_4->value();
    element[3] = ui->doubleSpinBox_5->value();
    element[4] = ui->doubleSpinBox_6->value();
    element[5] = ui->doubleSpinBox_7->value();
    element[6] = ui->doubleSpinBox_8->value()*M_PI/180;
    element[7] = ui->doubleSpinBox_9->value();

    element[0]=M2/M1;

    QString input2=ui->lineEdit_4->text();
    string data2 = input2.toUtf8().constData();
    std::ostringstream dat2NameStream(data2);
    dat2NameStream<<path<<"/"<<data2;
    std::string dat2Name = dat2NameStream.str();
    ifstream dat2(dat2Name.c_str());

    int lines=0;

    while(std::getline(dat2, line))
       ++ lines;

    dat2.clear();
    dat2.seekg(0, ios::beg);

    phim.resize(lines);
    mag.resize(lines);

    for (int i =0; i < lines; i++){
    dat2 >> eins >> zwei;
    istringstream istr(eins);
    istr >> phim[i];
    istringstream istr2(zwei);
    istr2 >> mag[i];
    //istringstream istr3(drei);
    //istr3 >> error[i];
    }
    dat2.close();

    string log = "logfile.dat";
    std::ostringstream logNameStream(log);
    logNameStream<<path<<"/"<<log;
    std::string logName = logNameStream.str();
    ofstream LogFile(logName.c_str());

    int nu = 8, eval=0;

    double y[nu+1], P[nu+1][nu], Z[nu], C[nu], S[nu], E[nu], e[nu+1][nu];

    int Gamma=2.0;
    int alpha=1.0;
    double beta=0.5;
    double btot=0.5;

    double step=0.05;

    for(int i =0; i<nu; i++){
        P[0][i]=element[i];
        //cout<<element[i]<<endl;
    }

    for(int i=0; i<nu+1; i++){
        for(int j=0; j<nu; j++){

            if((i>0)&(i==j+1)){
                e[i][j]=step;
            }
            else{
                e[i][j]=0;
            }
        }
    }

    for (int i=0; i<nu+1; i++){
        for (int j=0; j<nu; j++){

            P[i][j]=P[0][j]+e[i][j];
            element[j]=P[i][j];
            //cout<<element[j]<<" ";
        }
        eval++;
        qApp->processEvents(QEventLoop::AllEvents);
        y[i]=MainWindow::function();
        qApp->processEvents(QEventLoop::AllEvents);

        cout<<y[i]<<endl;
        cout<<endl;
        LogFile<<"Residuum of point "<<i<<": "<<y[i]<<endl;
    }

    int zaehler=40;
        int Ph, Pl, Psh;
    double yh, ysh, yl, ym, yi, ys, yt;

    //start main loop
        for (int t=0; t<zaehler; t++){

            LogFile<<endl;
            LogFile<<"Iteration: "<<t+1<<endl;
            qApp->processEvents(QEventLoop::AllEvents);

            cout<<endl;
            cout<<endl;
            cout<<"Iteration: "<<t+1<<"\t";

        //initialize next step
        ym=0;
        ys=0;
        for (int i=0; i<nu; i++){
        Z[i]=0;
        }

        //looking for highest value
        yh=y[0];
        Ph = 0;
        for (int j=0; j<nu+1; j++){
        if(y[j]>yh){
        yh = y[j];
        Ph = j;
        }}

        //looking for smallest value
        yl=yh;
        for (int j=0; j<nu+1; j++){
        if(y[j]<yl){
        yl=y[j];
        Pl = j;
        }}

        //looking for second highest value
        ysh=yl;
        for (int j=0; j<nu+1; j++){
        if((y[j]>ysh) & (y[j]<yh) & (y[j]>yl) & (j !=Pl)){
        ysh=y[j];
        Psh=j;
        }}

        yh=y[Ph];
        yl=y[Pl];
        ysh=y[Psh];

        //computing mean and sigma
        for (int i=0; i<nu+1; i++){
        ym+=y[i]/(nu+1);
        }
        for (int i=0; i<nu+1; i++){
        ys+=sqrt(pow((y[i]-ym),2));
        }
        ys=ys/(nu);

        cout<<"Simplex Mean: "<<ym<<" Simplex STD: "<<ys<<endl;

        LogFile<<"Mean: "<<ym<<"\tSTD: "<<ys<<endl;

        //compute centroid
        for (int j=0; j<nu; j++){
        for (int i=0; i<nu+1; i++){
        if (i!=Ph){
        Z[j]+=P[i][j]/nu;
        }}}

        //reflect highest value at centroid
        for (int i=0; i<nu; i++){
        C[i]=Z[i]+alpha*(Z[i]-P[Ph][i]);
                element[i]=C[i];
        }
        yi=MainWindow::function();
        eval++;

        if(yi<yl){
        for (int i=0; i<nu; i++){
        E[i]=Z[i]+Gamma*(C[i]-Z[i]);
                element[i]=E[i];
        }
        yt=MainWindow::function();
        eval++;
        if(yt<yl){
        for (int i=0; i<nu; i++){
        P[Ph][i]=E[i];
                element[i]=E[i];
        }
        y[Ph]=yt;//function(E);
        //eval++;
        }
        if (yt>=yl){
        for (int i=0; i<nu; i++){
        P[Ph][i]=C[i];
                element[i]=C[i];
        }
        eval++;
        y[Ph]=MainWindow::function();
        }}

        if(yi>=yl){
        if(yi<=ysh){
        for(int i=0; i<nu; i++){
        P[Ph][i]=C[i];
                element[i]=C[i];
        }
        eval++;
        y[Ph]=MainWindow::function();
        }
        if(yi>ysh){
        if(yi<=yh){
        for(int i=0; i<nu; i++){
        P[Ph][i]=C[i];
        }
        eval++;
        y[Ph]=MainWindow::function();
        yh=y[Ph];
        }
        for(int i=0; i<nu; i++){
        S[i]=Z[i]+beta*(P[Ph][i]-Z[i]);
        element[i]=S[i];
        }
        yt=MainWindow::function();
        eval++;
        if(yt>yh){
        for (int j=0; j<nu+1; j++){
        for (int i=0; i<nu; i++){
        P[j][i]=P[Pl][i]+btot*(P[j][i]-P[Pl][i]); //total contraction
        element[i]=P[j][i];
        }
        y[j]=MainWindow::function();
        eval++;
        }}

        if(yt<=yh){
        for(int i=0; i<nu; i++){
        P[Ph][i]=S[i];
        }
        eval++;
        y[Ph]=MainWindow::function();
        }}

        }

        }

}



double MainWindow::function()
{
    residu=0;

    for(int i =0; i<phim.size(); i++){

        q=element[0];
        R1=element[1];
        R2=element[2];
        L1=element[3];
        L2=element[4];
        R=element[5];
        inc=element[6];
        K=element[7];

        X=R*sin(2*M_PI*phim[i]);
        X1=-X/(1.0+(1.0/q));
        X2=X/(1.0+q);
        Y=R*cos(inc)*cos(2*M_PI*phim[i]);
        Y1=-Y/(1+(1/q));
        Y2=Y/(1+q);
        p=sqrt(pow((X2-X1),2)+pow((Y2-Y1),2));
        Z=R*sin(inc)*cos(2*M_PI*phim[i]);
        Z1=-Z/(1.0+(1.0/q));
        Z2=Z/(1.0+q);
        Theta1=2*acos((R2*R2-R1*R1-p*p)/(-2*R1*p));
        Theta2=2*acos((R1*R1-R2*R2-p*p)/(-2*R2*p));

        if(Theta1>M_PI){
            Theta1 = 2*M_PI-Theta1;
        }

        if(Theta2>M_PI){
            Theta2 = 2*M_PI-Theta2;
        }

        dA1=R1*R1*0.5*(Theta1-sin(Theta1));
        dA2=R2*R2*0.5*(Theta2-sin(Theta2));

        //cout<<Theta1<<"\t"<<Theta2<<"\t";

        if(Z1>Z2){
            if(p>R1+R2){
                A1 = R1*R1*M_PI;
                A2 = R2*R2*M_PI;
            }
            if((p<R1+R2) & (p>sqrt(R1*R1-R2*R2))){
                A1 = R1*R1*M_PI;
                A2 = R2*R2*M_PI-dA1-dA2;
            }
            if((p<sqrt(R1*R1-R2*R2)) & (p>R1-R2)){
                A1 = R1*R1*M_PI;
                A2 = dA2-dA1;
            }
            if(p<=R1-R2){
                A1 = R1*R1*M_PI;
                A2 = 0;
            }
        }
        if(Z1<Z2){
            if(p>R1+R2){
                A1 = R1*R1*M_PI;
                A2 = R2*R2*M_PI;
            }
            if((p<R1+R2) & (p>sqrt(R1*R1-R2*R2))){
                A1 = R1*R1*M_PI-dA1-dA2;
                A2 = R2*R2*M_PI;
            }
            if((p<sqrt(R1*R1-R2*R2)) & (p>R1-R2)){
                A1 = R1*R1*M_PI-R2*R2*M_PI+dA2-dA1;
                A2 = R2*R2*M_PI;
            }
            if(p<=R1-R2){
               A1 = R1*R1*M_PI-R2*R2*M_PI;
               A2 = R2*R2*M_PI;
            }
        }
        residu+=pow((mag[i]-K*(L1*A1/(R1*R1)+L2*A2/(R2*R2))/(4*M_PI)),2);
    }
    if(residu2==-1){
        residu2=residu;
    }
    else{
        if(residu<residu2){
            residu2=residu;
            cout<<q<<"\t"<<R1<<"\t"<<R2<<"\t"<<L1<<"\t"<<L2<<"\t"<<R<<"\t"<<inc*180/M_PI<<"\t"<<K<<"\t"<<residu<<endl;
            ui->doubleSpinBox_3->setValue(R1);
            ui->doubleSpinBox_4->setValue(R2);
            ui->doubleSpinBox_5->setValue(L1);
            ui->doubleSpinBox_6->setValue(L2);
            ui->doubleSpinBox_7->setValue(R);
            ui->doubleSpinBox_8->setValue(inc*180/M_PI);
            ui->doubleSpinBox_9->setValue(K);
            ui->doubleSpinBox_10->setValue(q);
        }
    }

    residu=residu/phim.size();

    return residu;
}

void MainWindow::on_doubleSpinBox_valueChanged()
{
    M1 = ui->doubleSpinBox->value();
    M2 = ui->doubleSpinBox_2->value();

    q=M2/M1;

    ui->doubleSpinBox_10->setValue(q);
}

void MainWindow::on_doubleSpinBox_2_valueChanged()
{
    M1 = ui->doubleSpinBox->value();
    M2 = ui->doubleSpinBox_2->value();

    q=M2/M1;

    ui->doubleSpinBox_10->setValue(q);
}
