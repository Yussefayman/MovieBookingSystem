
package moviebookingsystem;

import database.mysql;
import java.io.IOException;
import java.net.URL;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Random;
import java.util.ResourceBundle;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import javafx.scene.control.Label;


public class Final_winController implements Initializable {
     @FXML Label city;
     @FXML Label cinema;
     @FXML Label movie_name;
      @FXML Label time;
       @FXML Label date;
        @FXML Label chair;
         @FXML Label price;
          @FXML Label code;
           @FXML Label type;
           
    String city_= Submit_c_cController.city_;
    String cinema_= Submit_c_cController.cinema_;
    String movie_name_= FilmsController.movieName;
    String show_time = Date_time_submitController.time_;
   public static String show_date= Date_time_submitController.date_.toString();
   public static String show_type = Date_time_submitController.type_;
    String username= RegisterController.username_;
    String pass= RegisterController.pass;
    String u_n= LoginController.username_;
    String p = LoginController.pass;
    public static int cust_id;
     ArrayList<String>chairs = hall.chairs;
    public static String price_;
    public static int m_id;
    public static int cinema_id;
      static int counter=0;
    @FXML
    public void gen_code(ActionEvent event){
        counter++;
        if(counter==1){
            Random ran = new Random();
            int res =1000000 + ran.nextInt(999999);
            code.setText(Integer.toString(res));
        }
    }
   
    @Override
    public void initialize(URL url, ResourceBundle rb) {
     if(show_type.equals("2D")&&show_time.equals("15:00:00")){
         price_="30";
     }  
     else if(show_type.equals("2D")&&show_time.equals("18:00:00")){
          price_="40";
    }
     else if(show_type.equals("2D")&&show_time.equals("21:00:00")){
          price_="50";
    }
     else if(show_type.equals("3D")&&show_time.equals("00:00:00")){
          price_="70";
    }
     else{
         price_="100";
     }
     
     city.setText(city_);
     cinema.setText(cinema_);
     movie_name.setText(movie_name_);
     time.setText(show_time);
     date.setText(show_date);
     chair.setText(chairs.toString());
     price.setText(price_);
     type.setText(show_type);
    }    
     @FXML
    public void confirm_all(ActionEvent event) throws SQLException{
        Customer c= new Customer();
        cinema ci= new cinema();
        movie m=new movie();
        hall h= new hall();
        int h_id= h.hall_id(cinema_,city_);
         m_id=m.movie_id(movie_name_);
         cinema_id= ci.cinema_id(cinema_, city_);
        int show_id= h.show_id(show_time, show_type);
        /*
       if(isNullOrEmpty(username)&&isNullOrEmpty(pass)){
            cust_id= c.customer_id(u_n, p);
           for(int i=0;i<chairs.size();i++){
           mysql.runnonquery("insert into customer_hist (  customer_customer_id, chair_id, hall_id, movie_id, customer_hist_date , Ticket_Price,cinema_id,show_id,show_type) values ('"+cust_id+"','"+chairs.get(i)+"','"+h_id+"','"+m_id+"','"+show_date+"','"+price_+"','"+cinema_id+"','"+show_id+"','"+show_type+"')");
           }
       }
       else if(isNullOrEmpty(u_n)&&isNullOrEmpty(p)){
           int c_id= c.customer_id(username, pass);
           for(int i=0;i<chairs.size();i++){
           mysql.runnonquery("insert into customer_hist (  customer_customer_id, chair_id, hall_id, movie_id, customer_hist_date , Ticket_Price,cinema_id,show_id,show_type) values ('"+c_id+"','"+chairs.get(i)+"','"+h_id+"','"+m_id+"','"+show_date+"','"+price_+"','"+cinema_id+"','"+show_id+"','"+show_type+"')");
           }
       }*/
       
        Controls.callAlert("Info", "CONFIRM INFORMATION", "Your booking has been confirmed successfully");
    }
    @FXML 
    public void history(ActionEvent event) throws IOException{
        Controls control = new Controls();
           control.switchScenes("history");
    }
}
