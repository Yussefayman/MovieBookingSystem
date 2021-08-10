package moviebookingsystem;

import database.mysql;

import java.awt.image.BufferedImage;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URL;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ResourceBundle;
import java.util.logging.Level;
import java.util.logging.Logger;
import javafx.embed.swing.SwingFXUtils;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.fxml.FXMLLoader;
import javafx.fxml.Initializable;
import javafx.scene.Parent;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.stage.Stage;
import javax.imageio.ImageIO;




public class FilmsController implements Initializable  {
    
  @FXML Button p1 ;
  @FXML Button p2 ;
  @FXML Button p3 ;
  @FXML Button p4 ;
  @FXML Button p5 ;
  @FXML Button p6 ;
  @FXML ImageView im1;
  @FXML ImageView im2;
  @FXML ImageView im3;
  @FXML ImageView im4;
  @FXML ImageView im5;
  @FXML ImageView im6;
  @FXML Label l1;
  @FXML Label l2;
  @FXML Label l3;
  @FXML Label l4;
  @FXML Label l5;
  @FXML Label l6;
  @FXML Button showPoster ;
   
    public static String movieName;
    int mov_id_arr[];
    int movieSelectedID;
    int cinema_id = 0;
    String movie_poster[];
    String movies_name[];
  
     @FXML void selectPoster(ActionEvent e) throws IOException{
         
         Controls control = new Controls ();
         
          
        if(e.getSource()==p1){
            movieName=l1.getText();
            movieSelectedID=mov_id_arr[0];
            cinema_id= getCinemaId();
            FXMLLoader loader= new FXMLLoader();
            loader.setLocation(getClass().getResource("poster_info.fxml"));
            loader.load();
            Poster_infoController poster = loader.getController();
            poster.setMOVIEID(movieSelectedID);
             Parent p= loader.getRoot();
            Stage stage = new Stage ();
            stage.setScene(new Scene(p));
            stage.showAndWait();
           
        }
         if (e.getSource()== p2 ){
            movieName=l2.getText();
            movieSelectedID=mov_id_arr[1];
            cinema_id= getCinemaId();
            FXMLLoader loader= new FXMLLoader();
            loader.setLocation(getClass().getResource("poster_info.fxml"));
            loader.load();
            Poster_infoController poster = loader.getController();
            poster.setMOVIEID(movieSelectedID);
             Parent p= loader.getRoot();
            Stage stage = new Stage ();
            stage.setScene(new Scene(p));
            stage.showAndWait();
           
         }
         if (e.getSource()== p3 ){
            movieName=l3.getText();
           movieSelectedID=mov_id_arr[2];
            cinema_id= getCinemaId();
            FXMLLoader loader= new FXMLLoader();
            loader.setLocation(getClass().getResource("poster_info.fxml"));
            loader.load();
            Poster_infoController poster = loader.getController();
            poster.setMOVIEID(movieSelectedID);
             Parent p= loader.getRoot();
            Stage stage = new Stage ();
            stage.setScene(new Scene(p));
            stage.showAndWait();
         }
         if (e.getSource()== p4 ){
            movieName=l4.getText();
            movieSelectedID=mov_id_arr[3];
            cinema_id= getCinemaId();
            FXMLLoader loader= new FXMLLoader();
            loader.setLocation(getClass().getResource("poster_info.fxml"));
            loader.load();
            Poster_infoController poster = loader.getController();
            poster.setMOVIEID(movieSelectedID);
             Parent p= loader.getRoot();
            Stage stage = new Stage ();
            stage.setScene(new Scene(p));
            stage.showAndWait();
         }
         if (e.getSource()== p5 ){
            movieName=l5.getText();
            movieSelectedID=mov_id_arr[4];
            cinema_id= getCinemaId();
            FXMLLoader loader= new FXMLLoader();
            loader.setLocation(getClass().getResource("poster_info.fxml"));
            loader.load();
            Poster_infoController poster = loader.getController();
            poster.setMOVIEID(movieSelectedID);
             Parent p= loader.getRoot();
            Stage stage = new Stage ();
            stage.setScene(new Scene(p));
            stage.showAndWait();
         }
         if (e.getSource()== p6 ){
            movieName=l6.getText();
            movieSelectedID=mov_id_arr[5];
            cinema_id= getCinemaId();
            FXMLLoader loader= new FXMLLoader();
            loader.setLocation(getClass().getResource("poster_info.fxml"));
            loader.load();
            Poster_infoController poster = loader.getController();
            poster.setMOVIEID(movieSelectedID);
             Parent p= loader.getRoot();
            Stage stage = new Stage ();
            stage.setScene(new Scene(p));
            stage.showAndWait();
         }
            
     }
    
     @FXML
     public void showPoster(ActionEvent e) throws Exception{
         Button [] buttonArr = {p1,p2,p3,p4,p5,p6};
        ImageView [] imgs = {im1,im2 ,im3,im4,im5,im6};
        Label [] labelArr ={l1,l2,l3,l4,l5,l6};
        
          mov_id_arr = new int[7];
          movie_poster = new String[7];
          movies_name = new String [7];
          String city =l1.getText();
          String cinema = l2.getText();
          
          
          ResultSet rs = mysql.RunOnQuery("SELECT cinema_id FROM cinema WHERE cinema_name = '"+cinema+"' AND cinema_city = '"+city+"' ");
          while(rs.next()){
              cinema_id = rs.getInt("cinema_id");
          }
          ResultSet rs2 = mysql.RunOnQuery("SELECT movie_id FROM hall WHERE cinema_id= '"+cinema_id+"'");
          int idCount=0;
          while(rs2.next()){
              int movieID =rs2.getInt("movie_id");
              
              mov_id_arr[idCount]=movieID;
              System.out.print(mov_id_arr[idCount]);
              idCount++;
              
           
             
          }
          int counter =0 ;
          for(int i = 0  ; i < idCount ; i++){
              ResultSet rs3 = mysql.RunOnQuery("SELECT movie_name , movie_poster FROM movie Where movie_id = '"+mov_id_arr[i]+"' ");
              while(rs3.next()){
                  movie_poster[counter]=rs3.getString("movie_poster");
                  movies_name[counter]=rs3.getString("movie_name");
                  counter++;
              }
          }

        int count =0;
        for(int i = 0 ; i < counter ; i++){
            System.out.print(mov_id_arr[i]);
            String myPath ="C:\\Users\\Yousef\\Desktop\\pics\\";
            System.out.print(counter);
             myPath =myPath.concat(movie_poster[i]);
             System.out.print(myPath+"\n");
           FileInputStream in = new FileInputStream(myPath);
              BufferedImage img = ImageIO.read(in);
              Image image= SwingFXUtils.toFXImage(img, null);
              imgs[count].setImage(image);
              imgs[count].setVisible(true);
              buttonArr[count].setVisible(true);
              
              labelArr[count].setText(movies_name[count]); 
              labelArr[count].setVisible(true);
              count++;
        }
       
        
     }
       public int getCinemaId(){
         return  cinema_id;
     }
     public void setName(String city , String cinema){
         this.l1.setText(city);
         this.l2.setText(cinema);
     }
         
     

    @Override
    public void initialize(URL location, ResourceBundle resources) {
      


     
      } 
     

}
        
           

    
           
     
    

    


