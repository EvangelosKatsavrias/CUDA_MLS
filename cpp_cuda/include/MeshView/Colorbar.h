
#include<QApplication>

class Colorbar: public QWidget
{


public:

	Colorbar( QWidget* parent = 0 ) {

		this->setFixedHeight(30);

		QLinearGradient colorbar_gradient( 0, 0, this->width(), 0 );
		colorbar_gradient.setColorAt( 0, Qt::green );
		colorbar_gradient.setColorAt( 0.5, Qt::yellow );
		colorbar_gradient.setColorAt( 1, Qt::red );
		QPalette colorbar_palette;
		colorbar_palette.setBrush( QPalette::Background, colorbar_gradient );
		colorbar_palette.setBrush( QPalette::Text, Qt::white );
		this->setAutoFillBackground( true );
		this->setPalette( colorbar_palette );

/*		QHBoxLayout *layout = new QHBoxLayout(this);
		QLabel *minValue = new QLabel( QString( "-1" ), this );
		QLabel *zeroValue = new QLabel( QString( "0" ), this );
		QLabel *maxValue = new QLabel( QString( "1" ), this );

		QFont font = minValue->font();
		font.setPointSize( 10 ); 

		minValue->setFont(font); minValue->setFixedHeight( 20 );
		maxValue->setFont(font); maxValue->setFixedHeight( 20 );
		zeroValue->setFont(font);zeroValue->setFixedHeight( 20 );
*/
	//	layout->setContentsMargins( 0, 0, 0, 0 );
//		layout->setSizeConstraint( QLayout::SetMaximumSize );
//		QSpacerItem *space1 = new QSpacerItem( 10,100, QSizePolicy::Maximum );
//		QSpacerItem *space2 = new QSpacerItem( 10,100, QSizePolicy::Maximum );

//		layout->addStretch(0);
//		layout->addWidget(minValue);
//		layout->insertSpacerItem(1, space1);
//		layout->addWidget(zeroValue);
//		layout->insertSpacerItem(2, space2);
//		layout->addWidget(maxValue);

//		this->setLayout(layout);
	}


	void resizeEvent( QResizeEvent *event )
	{
		QWidget::resizeEvent( event );

		QLinearGradient colorbar_gradient( 0, 0, this->width(), 0 );
		colorbar_gradient.setColorAt( 0, Qt::green );
		colorbar_gradient.setColorAt( 0.5, Qt::yellow );
		colorbar_gradient.setColorAt( 1, Qt::red );
		QPalette colorbar_palette;
		colorbar_palette.setBrush( QPalette::Background, colorbar_gradient );
		colorbar_palette.setBrush( QPalette::Text, Qt::white );
		this->setPalette( colorbar_palette );

	}
};



class QualityColorbar: public QWidget
{


public:

	QualityColorbar( QWidget* parent = 0 ) {

		this->setFixedHeight(30);

		QLinearGradient colorbar_gradient( 0, 0, this->width(), 0 );
		colorbar_gradient.setColorAt( 0, Qt::green );
		colorbar_gradient.setColorAt( 0.5, Qt::yellow );
		colorbar_gradient.setColorAt( 1, Qt::red );
		QPalette colorbar_palette;
		colorbar_palette.setBrush( QPalette::Background, colorbar_gradient );
		colorbar_palette.setBrush( QPalette::Text, Qt::white );
		this->setAutoFillBackground( true );
		this->setPalette( colorbar_palette );

/*		QHBoxLayout *layout = new QHBoxLayout(this);
		QLabel *minValue = new QLabel( QString( "-1" ), this );
		QLabel *zeroValue = new QLabel( QString( "0" ), this );
		QLabel *maxValue = new QLabel( QString( "1" ), this );

		QFont font = minValue->font();
		font.setPointSize( 10 ); 

		minValue->setFont(font); minValue->setFixedHeight( 20 );
		maxValue->setFont(font); maxValue->setFixedHeight( 20 );
		zeroValue->setFont(font);zeroValue->setFixedHeight( 20 );
*/
	//	layout->setContentsMargins( 0, 0, 0, 0 );
//		layout->setSizeConstraint( QLayout::SetMaximumSize );
//		QSpacerItem *space1 = new QSpacerItem( 10,100, QSizePolicy::Maximum );
//		QSpacerItem *space2 = new QSpacerItem( 10,100, QSizePolicy::Maximum );

//		layout->addStretch(0);
//		layout->addWidget(minValue);
//		layout->insertSpacerItem(1, space1);
//		layout->addWidget(zeroValue);
//		layout->insertSpacerItem(2, space2);
//		layout->addWidget(maxValue);

//		this->setLayout(layout);
	}


	void resizeEvent( QResizeEvent *event )
	{
		QWidget::resizeEvent( event );

		QLinearGradient colorbar_gradient( 0, 0, this->width(), 0 );
		colorbar_gradient.setColorAt( 0, Qt::green );
		colorbar_gradient.setColorAt( 0.375, Qt::cyan );
		colorbar_gradient.setColorAt( 0.65, Qt::yellow );
		colorbar_gradient.setColorAt( 0.875, Qt::magenta );
		colorbar_gradient.setColorAt( 0.97, Qt::red );
		colorbar_gradient.setColorAt( 1, Qt::darkRed );
		QPalette colorbar_palette;
		colorbar_palette.setBrush( QPalette::Background, colorbar_gradient );
		colorbar_palette.setBrush( QPalette::Text, Qt::white );
		this->setPalette( colorbar_palette );

	}
};


