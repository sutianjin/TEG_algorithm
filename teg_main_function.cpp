
//本文件主要是说明TEG上位机应用软件处理采样数据，
//主要包括了数据点的均值滤波、曲线轮廓提取、主要检测参数提取

//存储滤波后的数据
QList<QPointF> List_Filter_Data；

//存储凝血曲线轮廓的上支曲线
QList<QPointF> List_Max_Curve_Data；

//存储凝血曲线轮廓的下支曲线
QList<QPointF> List_Max_Curve_Data；

/*
@数据点采样均值滤波处理：下位机按照5ms发送一个采样点
@入参：orign_data,串口采样的原始数据
@      size,采样数据大小
@出参：无
@返回值：无
*/
void Smooth_Origin_Data(double* orign_data, unsigned size)
{
	//如果采样的原始数据大小不是36的整数倍，则继续采样，因此在传之前最好先判断是否为36的整数倍
	if(size%36!=0) return;
	
	//采用均值滤波，N值选择为36，因为采用的一个小周期正好是36个点
	for(unsigned i = 0; i < size; i+=36)
	{
		int sum = 0;
		int ave = 0;
		unsigned int j;
		for(j = i; j< i+36; j++){

			sum += data[j];
		}
		ave = sum/36;
		
		//主要是存储滤波后的数据点，选择中值作为滤波后的数据点
		QPointF s(double(j+18+List_Filter_Data.size())*0.005,ave);	
		List_Filter_Data.append(s);
	}	
}

/*
@凝血曲线的轮廓提取
@入参：List_Filter_Data，滤波后的数据
@出参：无
@返回值：无
*/
void ExtractionContour(QList<QPointF> List_Filter_Data)
{
	unsigned int size = List_Filter_Data.size();
	//由于TEG旋转一个周期为10s左右，滤波后的一个正弦波周期也为10s， 10000/36 = 55.5，故选择55个点为一个轮廓点的周期
	for(unsigned int i = 0; i < size-55,i=i+55)
	{
		int xmax = 0;
		int xmin = 0;
		int ymax = 0;
		int ymin = 5000;

		for(unsigned int j = i;j<i-55;j++){
			if(List_Filter_Data.at(j).ry()> ymax){
				ymax = List_Filter_Data.at(j).ry();
				xmax = j-1;
			}
			if(List_Filter_Data.at(j).ry()< ymin){
				ymin = List_Filter_Data.at(j).ry();
				xmin = j-1;
			}
		
		QPointF max_point(xmax,ymax);
		QPointF min_point(xmin,ymin);
		List_Max_Curve_Data.append(max_point)；
		List_Min_Curve_Data.append(min_point)；
		
	}	
}
 
 
/*
@获取主要的凝血检测参数
@入参：List_Max_Curve_Data，曲线上支轮廓
@      List_Min_Curve_Data，曲线下值轮廓
@出参：R_Point,R值坐标点
@	   K_Point,K值坐标点
@      Angle_Point,angle角度值坐标点
@      MA_Point,MA值坐标点
@返回值：无
*/

void GetMainParameter(QList<QPointF> List_Max_Curve_Data,QList<QPointF> List_Min_Curve_Data，QPointF &R_Point,QPointF &K_Point,QPointF &Angle_Point,QPointF &MA_Point)
{
	MA_Point = QPointF(0,5000);
	unsigned int size = List_Max_Curve_Data.size();
	for(unsigned int i = 0; i < size; ++i)
	{
		//R值表示了TEG振荡幅值超过20mm为R值点。
		if(List_Max_Curve_Data[i]>=(List_Max_Curve_Data[i]+List_Min_Curve_Data[i])/2 + 20){
			//K值表示TEG振荡幅值超过500的点
			if(List_Max_Curve_Data[i]<=(List_Max_Curve_Data[i]+List_Min_Curve_Data[i])/2 + 500){
				K_Point = List_Max_Curve_Data[i];
			}
			
			//设置MA值,
			if(List_Min_Curve_Data[i].ry() < MA_Point.ry()){
				MA_Point = List_Min_Curve_Data[i];
			}			
			
			//求angle角坐标
			lineParaPoint = linefunc(RvaluePoint,List_Max_Curve_Data[i]);
			if((List_Max_Curve_Data[i-2].y() < getlineY(lineParaPoint,List_Max_Curve_Data[i-2].x()))
				&&((List_Max_Curve_Data[i-3].y() < getlineY(lineParaPoint,List_Max_Curve_Data[i-3].x())))){
				AnglevaluePoint = List_Max_Curve_Data[i-1];
			}
			
		}else{
			//如果振荡幅值不超过20mm，则将其规划为一条平稳线，另外最大轮廓值也可以表示振荡的实时值。
			List_Max_Curve_Data.append(QPointF(List_Max_Curve_Data[i].rx(),List_Max_Curve_Data[i]+List_Min_Curve_Data[i])/2);
			List_Max_Curve_Data.append(QPointF(List_Min_Curve_Data[i].rx(),List_Min_Curve_Data[i]+List_Min_Curve_Data[i])/2);
			R_Point = List_Max_Curve_Data[i];
		}
	}			
}


/*
@获取两点的斜率k和截距b ；y = kx+b
@入参：point1,第一个点
	   point2,第二个点
@出参：无
@返回值：QPointF(k,b)
*/
QPointF linefunc(QPointF point1, QPointF point2)
{
    QPointF linepara;
    linepara.setX((point1.ry()-point2.ry())/(point1.rx()-point2.rx()));
    linepara.setY(point1.ry()-(linepara.rx()*point1.rx()));
    return linepara;
}

/*
@获取曲线的Y轴值
@入参：lineSlope，表示曲线的斜率和截距
	   xValue,X轴的值
@出参：无
@返回值：Y轴的值
*/
double getlineY(QPointF lineSlope, double xValue)
{
    double  yValue;
    yValue = (lineSlope.rx()*xValue)+lineSlope.ry();
    return yValue;
}


/*
@设置R值
@入参： R_Point R值坐标点
@出参： strRvalue，R值的字符串值
返回值：无
*/
void SetRvalue(QPointF R_Point,QString &strRvalue)
{

    double value = R_Point.rx();
    QString strRvalue = QString::number(value/60,'f',2);
}


/*
@设置K值
@入参： K_Point K值坐标点
@出参： strKvalue，K值的字符串值
返回值：无
*/
void SetRvalue(QPointF K_Point,QString &strKvalue)
{
    double value = K_Point.rx();
    QString strKvalue = QString::number(value/60,'f',2);
}


/*
@设置MA值
@入参： MA_Point MA值坐标点
        baselineValue 基线的y轴值
@出参： strMAvalue,MA值的字符串值
返回值：无
*/
void SetRvalue(QPointF MA_Point,double baselineValue,QString &strMAvalue)
{
	//MA值是曲线轮廓下半支与基线之差*0.047
    double value = MA_Point.ry();
    QString strMAvalue = QString::number((baselineValue-value)*0.047,'f',2);
}

/*
@设置ANGLE值
@入参： ANGLE_Point ANGLE值坐标点
        lineParaPoint 曲线的斜率和截距点
@出参： strAnglevalue,MA值的字符串值
返回值：无
*/
void SetAnglevalue(QPointF AnglevaluePoint, QPointF lineParaPoint,,QString &strAnglevalue)
{
	//angle角度值是切点与R值点的斜率。k=((yangle-yRvalue)*0.047)/((xangle-xRvalue)*4/60)
    double tanvalue = ((AnglevaluePoint.ry()-preRValuePoint[channel].ry())*0.047)/((AnglevaluePoint.rx()-preRValuePoint[channel].rx())*4/60);
    double arctanvalue = (180*qAtan(tanvalue))/M_PI;
    QString strAnglevalue = QString::number(arctanvalue,'f',2);
}














