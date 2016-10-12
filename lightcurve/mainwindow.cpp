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

using namespace std;

double M1, M2, R1, R2, L1, L2, R, q, X, Y, X1, X2, Y1, Y2, Z, Z1, Z2, inc, A1, A2, dA1, K, dA2, Phi, Theta1, Theta2, p;
QString QPath;
string path, line, eins;

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

    ui->plot->axisRect()->setupFullAxesBox(true);
    ui->plot->yAxis->setRange(0,1.2);
    ui->plot->xAxis->setRange(0,1.5);

    ui->lineEdit->setText("phases.dat");
    ui->lineEdit_2->setText("light.dat");
    ui->lineEdit_3->setText("/home/daniels/LightCurve");


}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_pushButton_2_clicked()
{
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

        I[i] = (L1*A1/(R1*R1)+L2*A2/(R2*R2))/(4*M_PI);

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
    ui->plot->addGraph();
    ui->plot->graph(0)->setData(phi, I);
    ui->plot->addGraph();
    ui->plot->graph(1)->setData(phi, I1);
    ui->plot->addGraph();
    ui->plot->graph(2)->setData(phi, I2);
    ui->plot->addGraph();
    ui->plot->graph(3)->setData(phi, Ratio);
    ui->plot->graph(0)->setPen(pen);
    ui->plot->graph(1)->setPen(pen2);
    ui->plot->graph(2)->setPen(pen3);
    ui->plot->graph(3)->setPen(pen4);

    ui->plot->replot();


}
