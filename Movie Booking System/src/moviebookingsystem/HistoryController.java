package moviebookingsystem;
import database.mysql;
import java.net.URL;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ResourceBundle;
import java.util.logging.Level;
import java.util.logging.Logger;
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import javafx.scene.control.Label;


public class HistoryController implements Initializable {
    @FXML Label chairs;
    @FXML Label movie;
    @FXML Label hall;
    @FXML Label date;
    @FXML Label price;
    @FXML Label cinema;
    @FXML Label type;
    int id= Final_winController.cust_id;
    int movie_id= Final_winController.m_id;
   
    @Override
    public void initialize(URL url, ResourceBundle rb) {
      int mov_id = 0,cinema_ = 0;
        try {
            ResultSet rs= mysql.RunOnQuery("select * from customer_hist where customer_customer_id= '"+id+"'");
             System.out.println(id);
            System.out.println(rs.next());
       while(rs.next()){
           String chair_id= rs.getString("chair_id");
            String hall_id= rs.getString("hall_id");
              mov_id = rs.getInt("movie_id");
            String date_= rs.getString("customer_hist_date");
            String price_ =rs.getString("Ticket_Price");
            cinema_= rs.getInt("cinema_id");
            String show_type=rs.getString("show_type");
            System.out.println(movie_id+" "+mov_id);
            if(movie_id==mov_id){
            chairs.setText(chair_id);
            hall.setText(hall_id);
            date.setText(date_);
            price.setText(price_);
            type.setText(show_type);
            }
       }
         ResultSet rs2= mysql.RunOnQuery("select * from movie");
       while(rs2.next()){  
           String movie_n=rs2.getString("movie_name");
           int mo_id= rs2.getInt("movie_id");
           if(mo_id==mov_id&&mo_id==movie_id)
           {
               movie.setText(movie_n);
               
               
           }
       }
        ResultSet rs3= mysql.RunOnQuery("select * from cinema");
       while(rs3.next()){
           int cinema_id= rs3.getInt("cinema_id");
           String cinema_n=rs3.getString("cinema_name");
          
           if( cinema_id==cinema_)
           {
              cinema.setText(cinema_n);
               
           }
       }
       
       
}       catch (SQLException ex) {
            Logger.getLogger(HistoryController.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}