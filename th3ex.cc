void th3ex(){

	TH3C *volume = new TH3C("Volume","Neutron production volume", 10,-50.,50.,10,-50.,50.,10,-50.,50.);
	TCanvas *canvas = new TCanvas("Canvas","Neutron production volume",600,600);

	float x,y,z;
	for(x=-50;x<50;x+=1){
		for(y=-50;y<50;y+=1){
			for(z=-50;z<50;z+=1){
			if((x*x+y*y+.3*z*z)<500)
			volume->Fill(x,y,z,1);
			}
		}
	}
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
			volume->Fill(-50.,-50.,0.,1);
	
	gStyle->SetCanvasPreferGL(true);
	canvas->Draw();
	volume->Draw("glcolz");
}
